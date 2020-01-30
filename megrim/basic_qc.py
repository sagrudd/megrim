#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:11:57 2019

@author: srudd
"""

import logging
import pandas as pd
import numpy as np
import math
import sys
import os
from math import log10
import matplotlib as mpl
from tqdm import tqdm
from scipy import stats
from megrim.environment import Flounder
from bisect import bisect_left
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from time import time
from megrim.infographic_plots import InfographicPlot, InfographicNode
from bokeh.models import LinearColorMapper, BasicTicker, ColorBar, Label, \
    LabelSet, NumeralTickFormatter, Span, ColumnDataSource, Range1d
from bokeh.palettes import (Blues9)
from bokeh.plotting import figure
import functools

# import palettable.colorbrewer.sequential


class Timer():
    """
    https://stackoverflow.com/questions/3620943/measuring-elapsed-time-with-the-time-module/46544199
    """

    def __init__(self, message):
        self.message = message

    def __enter__(self):
        self.start = time()
        return None

    def __exit__(self, type, value, traceback):
        e = int(time() - self.start)
        m = '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60)
        print(self.message.format(m))


class SequenceSummaryHandler(Flounder):

    def __init__(self, target_file=None, target_data=None, fcid=None):
        Flounder.__init__(self)
        self.target_file = target_file
        self.fcid = fcid
        if target_file is not None:
            self._import_data()
        elif target_data is not None:
            self.seq_sum = target_data
            self.seq_sum_head = target_data.head()


    def _load_seqsum(self, file=None):
        if file is None:
            file = self.target_file
        extension = os.path.splitext(file)[1].lower()
        self.seq_sum = None
        blocksize = 64000000
        compression = None
        if (extension == ".bz2"):
            logging.warning(
                "reading a bzip2 file - this has performance implications")
            blocksize = None
            compression = "bz2"
        elif (extension in [".gzip", ".gz"]):
            logging.warning(
                "reading a gzip file - this has performance implications")
            blocksize = None
            compression = "gzip"
        self.seq_sum = dd.read_csv(
            file,
            delimiter='\t',
            blocksize=blocksize,
            compression=compression
        )
        # slice out the head entries for further functionality
        self.seq_sum_head = self.seq_sum.head()


    def _import_data(self):
        """
        method parses the provided Sequencing_summary.txt file and loads 
        the reduced contents into memory - a DASK dataframe is used to
        allow for the scaling to large PromethION runs (or equivalent)

        There is a possibility of duplicate header lines; these will btreak
        the dask import
        
        grep can be used to identify the line numbers - this can take a while
        on a promethion flowcell

        Returns
        -------
        None.

        """

        try:
            self._load_seqsum()
        except ValueError as verr:
            logging.debug("ValueError error: {0}".format(verr))
            logging.warning(
                "Dask ValueError during import - filtering duplicate headers")

            extension = os.path.splitext(self.target_file)[1].lower()
            compression = None
            blocksize = 6400000
            if (extension == ".bz2"):
                compression = "bz2"
            elif (extension in [".gzip", ".gz"]):
                compression = "gzip"
            
            data = dd.read_csv(self.target_file, delimiter="\t", 
                               compression=compression, blocksize=None, 
                               dtype="object")
            data = data[~(data["filename"] == 'filename')].compute()
            type_info = {"channel": "int64",
                         'start_time': "float64",
                         'duration': "float64",
                         'num_events': "int64",
                         'sequence_length_template': "int64",
                         'mean_qscore_template': "float64",
                         }
            # older versions of albacore used True | recent versions use TRUE
            data.passes_filtering = data.passes_filtering.isin(
                ["True", "TRUE"])
            
            for key in type_info.keys():
                if key in data.columns:
                    data = data.astype({key: type_info.get(key)})
            self.seq_sum = dd.from_pandas(data, chunksize=blocksize)
            self.seq_sum_head = self.seq_sum.head()

        except Exception as e:
            logging.error("ValueError error: {0}".format(e))
            logging.error("ERROR - this is an unexpected edge case ...")
            sys.exit(0)
        # start excluding dask columns that are not of core interest
        keep = ['read_id', 'channel', 'start_time', 'duration', 'num_events',
                'sequence_length_template', 'mean_qscore_template',
                'passes_filtering', 'barcode_arrangement',
                ]
        for col in self.seq_sum.columns:
            if not col in keep:
                logging.debug("dropping %s" % col)
                self.seq_sum = self.seq_sum.drop(col, axis=1)
        pbar = ProgressBar()
        pbar.register()
        with Timer("Elapsed time to extract sequence data {}"):
            self.seq_sum = self.seq_sum.compute()
            self.seq_sum = self.seq_sum.set_index('read_id')
        pbar.unregister()

    @functools.lru_cache()
    def import_barcodes_file(self, barcode_file):
        logging.info("importing barcode file {}".format(barcode_file))

        extension = os.path.splitext(barcode_file)[1].lower()
        compression = None
        blocksize = 6400000
        if extension == ".bz2":
            compression = "bz2"
            blocksize = None
        elif extension in [".gzip", ".gz"]:
            compression = "gzip"
            blocksize = None

        pbar = ProgressBar()
        pbar.register()
        data = dd.read_csv(barcode_file, delimiter="\t", 
                           compression=compression, blocksize=blocksize, 
                           dtype="object")
        data = data.iloc[:,[0,1]].compute()
        pbar.unregister()
        data = data.set_index('read_id')
        return data
    
    def merge_barcodes_file(self, barcode_file):
        bc_data = self.import_barcodes_file(barcode_file)
        self.seq_sum = pd.concat([self.seq_sum, bc_data], axis=1, join='inner')
        

    def get_flowcell_id(self):
        """
        This method pulls the flowcell id from the fastq filename of the 
        original sequence ... this code may need to be maintained ...

        Returns
        -------
        fastq_id : TYPE
            DESCRIPTION.

        """

        columns = list(self.seq_sum_head.columns)
        target = None
        if "filename_fastq" in columns:
            target = "filename_fastq"
        elif "filename" in columns:
            target = "filename"
        else:
            logging.warning("there is not a suitable filename column")
            if self.fcid is not None:
                return self.fcid
            else:
                return "undefined"

        fastq_id = str(self.seq_sum_head[target].iloc[0]).split("_")
        if len(fastq_id) == 3:
            fastq_id = fastq_id[0]
        elif len(fastq_id) == 15:
            fastq_id = fastq_id[2]
        else:
            logging.warning(
                "ERROR - unable to parse flowcell_id from seq_summary ... "+
                str(fastq_id))
            fastq_id = "undefined"
        logging.debug("flowcell read as [%s]" % fastq_id)
        # print(("_", fastq_id))
        return fastq_id

    def get_read_count(self):
        return len(self.seq_sum.index)

    def executive_summary(self, **kwargs):
        """
        Method prepares a three panel infographic plot that summarises the
        flowcell, reads and bases sequenced within the associated run

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        (plot_width, plot_dpi) = self.handle_kwargs(
            ["plot_width", "plot_dpi"], **kwargs)
        read_count = self.get_read_count()
        total_bases = self.seq_sum['sequence_length_template'].sum()
        
        flowcell_node = InfographicNode(legend="flowcell",
                                        value=self.get_flowcell_id(),
                                        graphic='fingerprint')

        readcount_node = InfographicNode(legend="Reads produced",
                                         value="{:,}".format(read_count),
                                         graphic='filter')

        gb_seq_val = InfographicNode(legend="Gigabases called",
                                     value="{:.2f}".format(total_bases / 1e9),
                                     graphic='flag-checkered')
        infographic_data = [flowcell_node, readcount_node, gb_seq_val]
        ip = InfographicPlot(infographic_data, rows=1, columns=3)
        return ip.plot_infographic(plot_width, plot_dpi)

    def get_absolute_runtime(self):
        max_time = self.seq_sum['start_time'].max()
        return max_time

    def get_runtime(self, units="hours", rounded=True):
        """
        Accessory method to return the length of time that elapsed during 
        this flowcell run. If rounded (default) this will be rounded *up*
        to the nearest canonical timing ...

        Parameters
        ----------
        units : TYPE, optional
            DESCRIPTION. The default is "hours"
        rounded : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        runtime : TYPE
            DESCRIPTION.

        """
        scale = {"minutes": 60,
                 "hours": 60 * 60,
                 "days": 60 * 60 * 24}
        canonical_runs = [1, 2, 4, 6, 8, 12, 18, 24, 36, 48, 60, 72]

        if not units in scale.keys():
            raise Exception("{%s} is not a qualified runtime unit ...")

        runtime = math.ceil(self.get_absolute_runtime() / scale["hours"])
        if rounded:
            def take_closest(values, val):
                pos = bisect_left(values, val)
                if pos == len(values):
                    return values[-1]
                else:
                    return values[pos]

            runtime = take_closest(canonical_runs, runtime)
        return (runtime * scale["hours"]) / scale[units]

    def plot_passed_gauge(self, **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], 
            **kwargs)
        
        read_count = len(self.seq_sum)
        passed_read_count = self.seq_sum.passes_filtering.sum()
        perc_val = passed_read_count / read_count * 100

        p = figure(plot_width=plot_width, plot_height=plot_height, 
                   x_range=(0.25, 1.75), y_range=(0.7, 1.5), tools=plot_tools)

        start_val = 0
        middle_val = (math.pi / 100) * (100 - perc_val)
        end_val = math.pi

        p.annular_wedge(
            x=[1], y=[1], inner_radius=0.2, outer_radius=0.5, 
            start_angle=middle_val, end_angle=end_val, color="green", 
            alpha=0.6)

        p.annular_wedge(
            x=[1], y=[1], inner_radius=0.2, outer_radius=0.5,
            start_angle=start_val, end_angle=middle_val, color="orange", 
            alpha=0.6)

        label = Label(
            x=1, y=1, text="{:.1f}%".format(perc_val), x_units='data', 
            y_units='data', text_align='center', text_font_style='bold', 
            text_font_size='1.5em')
        legend = Label(x=1, y=0.9,
                       text="Percentage of reads passing QC filter",
                       x_units='data', y_units='data',
                       text_align='center', text_font_size='1.9em')
        p.add_layout(label)
        p.add_layout(legend)
        p.axis.visible = False
        p.xgrid.visible = False
        p.ygrid.visible = False
        return self.handle_output(p, plot_type, prefix="passed_gauge")


    def plot_channel_activity(self, **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        channel_map = SequencingSummaryGetChannelMap(self.seq_sum)
        # layout = channel_map.get_platform_map()
        layout = channel_map.get_platform_density()
        layout = layout.fillna(0)
        layout['row'] = layout['row'].astype(str)
        layout['column'] = layout['column'].astype(str)
        layout['count'] = layout['count'].astype(int)

        logging.debug(layout)

        colors = Blues9[::-1]
        mapper = LinearColorMapper(
            palette=colors, low=layout['count'].min(), 
            high=layout['count'].max())

        rows = list(layout.row.unique())
        columns = list(layout.column.unique())

        p = figure(
            title="channel activity plot", x_range=columns, y_range=rows,
            x_axis_location="above", plot_width=plot_width, 
            plot_height=plot_height, tools=plot_tools, 
            toolbar_location='below')

        p.axis.visible = False
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "5pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = math.pi / 3
        #p.title.text_font_size = '18pt'

        p.rect(x="column", y="row", width=1, height=1,
               source=layout,
               fill_color={'field': 'count', 'transform': mapper},
               line_color=None)

        color_bar = ColorBar(
            color_mapper=mapper, major_label_text_font_size="10pt",
            ticker=BasicTicker(desired_num_ticks=len(colors)),
            title="#reads", label_standoff=6, border_line_color=None, 
            location=(0, 0))
        p.add_layout(color_bar, 'right')
        return self.handle_output(p, plot_type)


    @functools.lru_cache()
    def get_passed_reads(self, field=None, passed=True):
        if passed is None:
            return self.seq_sum
        else:
            if field is None:
                return self.seq_sum.loc[
                    self.seq_sum['passes_filtering']==passed, ]
            else:
                return self.seq_sum.loc[
                    self.seq_sum['passes_filtering']==passed, field]


    @functools.lru_cache()
    def get_length_character(self, key, field="sequence_length_template", 
                             passed=None):
        """
        This is an accessory method to help get a sequence length threshold
        based on parameters that may be more interesting thanj just the
        canonical longest read or read mean length

        Parameters
        ----------
        key : TYPE
            The feature to report - this can include keywords (max), (mean),
            (qmean), (N<float>), (n<float>) or (Q<float). 
        field : TYPE, optional
            DESCRIPTION. The default is "sequence_length_template".
        passed : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        int
            a sequence length that may be used for downstream processing and
            visualisation.

        """
        if isinstance(key, int):
            return key
        elif key is None:
            return self.get_passed_reads(field=field, passed=passed).max()
        elif key == "max":
            return self.get_passed_reads(field=field, passed=passed).max()
        elif key == "mean":
            return self.get_passed_reads(field=field, passed=passed).mean()
        elif key == "min":
            return self.get_passed_reads(field=field, passed=passed).min()
        elif key == "qmean":
            series = self.get_passed_reads(field=field, passed=passed)
            return -10 * log10((10 ** (series / -10)).mean())
        # handle N.values
        elif key[0]=="N":
            val = float(key[1:])
            ldata = self.get_passed_reads(
                field=field, passed=passed).sort_values(
                    ascending=False).reset_index(drop=True)
            n_sum = int(ldata.sum())
            n_targ = (n_sum / 100) * val
            accumulated = ldata.cumsum()
            aindex = accumulated.loc[(accumulated >= n_targ)].index[0]
            return ldata[aindex]
        elif key[0]=="n":
            val = int(key[1:])
            ldata = self.get_passed_reads(
                field=field, passed=passed).sort_values(
                    ascending=False).reset_index(drop=True)
            return ldata[val]
        elif key[0]=="Q":
            val = float(key[1:])
            return self.get_passed_reads(
                field=field, passed=passed).quantile(val)
        logging.error("{} is not an understood transformation")
        return None


    def library_characteristics_infographic(self, **kwargs):

        (plot_width, plot_dpi) = self.handle_kwargs(
            ["plot_width", "plot_dpi"], **kwargs)
        
        longest_read = self.get_length_character(key="max", passed=True)
        mean_read_length = self.get_length_character(key="mean", passed=True)
        read_n50_length = self.get_length_character(key="N50", passed=True)

        passed_mean_q = self.get_length_character(
            key="qmean", passed=True, field="mean_qscore_template")
        failed_mean_q = self.get_length_character(
            key="qmean", passed=False, field="mean_qscore_template")

        mean_read_length_node = InfographicNode(
            legend="Mean Read Length", 
            value="{:.2f}".format(mean_read_length),
            graphic='map-signs')
        read_n50_length_node = InfographicNode(
            legend="N50",
            value="{:,}".format(read_n50_length),
            graphic='bullseye')
        passed_mean_q_node = InfographicNode(
            legend="Mean Read Quality",
            value="{:.2f}".format(passed_mean_q),
            graphic='award')
        failed_mean_q_node = InfographicNode(
            legend="Mean Failed QV",
            value="{:.2f}".format(failed_mean_q),
            graphic='bug')
        longest_read_node = InfographicNode(
            legend="Longest Read",
            value="{:,}".format(longest_read),
            graphic='sort')
        infographic_data = [mean_read_length_node, 
                            read_n50_length_node,
                            passed_mean_q_node, 
                            failed_mean_q_node,
                            longest_read_node]
        ip = InfographicPlot(infographic_data, rows=1, columns=5)
        return ip.plot_infographic(plot_width, plot_dpi)



    @functools.lru_cache()
    def extract_size_stratified_data(self, longest_read, bins):
        
        l_seq_sum = self.seq_sum[[
            "sequence_length_template", "passes_filtering"]]
        
        longest_read = self.get_length_character(key=longest_read, passed=True)
        logging.debug("longest_read == {}".format(longest_read))
        
        longest_read = int(longest_read + 1)
        boundaries = np.linspace(
            0, longest_read, num=bins, endpoint=True, retstep=False)
        assignments = np.digitize(
            l_seq_sum.sequence_length_template, boundaries)
        
        l_seq_sum = l_seq_sum.reindex(
            columns=l_seq_sum.columns.tolist()+
            ["batch", "pass_bases", "fail_bases", "pass_reads", "fail_reads"])
        
        l_seq_sum.iloc[:,[2]] = assignments
        l_seq_sum.iloc[:,[3,4,5,6]] = 0

        # There is a logic issue - at the time of writing it was assumed that
        # sequence collections would be native and contain a mix of pass
        # and fail sequences ...

        l_seq_sum.loc[~l_seq_sum.passes_filtering, ["fail_reads"]] = 1
        l_seq_sum.loc[l_seq_sum.passes_filtering, ["pass_reads"]] = 1
        l_seq_sum.loc[
            l_seq_sum.passes_filtering, ["pass_bases"]
            ] = l_seq_sum.loc[
                l_seq_sum.passes_filtering, ["sequence_length_template"]
                ].iloc[:,0]
        l_seq_sum.loc[~l_seq_sum.passes_filtering, ["fail_bases"]] = l_seq_sum.loc[~l_seq_sum.passes_filtering, ["sequence_length_template"]].iloc[:,0]
        
        l_seq_res = l_seq_sum.groupby(["batch"]).agg(
            {"batch": "first",  "fail_reads": np.sum, "pass_reads": np.sum,
             "pass_bases": np.sum, "fail_bases": np.sum })
    
        # handle the possibility of missing values ... e.g. Dolores
        l_seq_res = l_seq_res.reindex(pd.Index(pd.Series(boundaries).index, name="hh")).reset_index()
        l_seq_res = l_seq_res.drop("hh", axis=1)
    
        # shape the data for quad presentation ...
        l_seq_res = l_seq_res.reindex(columns=l_seq_res.columns.tolist()+["left", "right"])
        l_seq_res.iloc[1:,[5]] = boundaries[:-1, np.newaxis]
        l_seq_res.iloc[1:,[6]] = boundaries[1:, np.newaxis]
        return (longest_read,
                boundaries, 
                l_seq_res.fillna(0), 
                self.get_length_character(key="mean", passed=True), 
                self.get_length_character(key="N50", passed=True))


    def plot_sequence_length(self, normalised=True,
                             include_failed=True, bins=30,
                             annotate_mean=True, annotate_n50=True, 
                             longest_read="Q0.995", **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        longest_read, boundaries, l_seq_res, mean_val, N = self.extract_size_stratified_data(longest_read=longest_read, bins=bins)
        if include_failed:
            if normalised:
                failed = pd.DataFrame({"bases":l_seq_res["fail_bases"].tolist(), "bases_line": 0, "classification": "failed", "colour": "#A6CEE3", "left": l_seq_res["left"].tolist(), "right": l_seq_res["right"].tolist()})
                passed = pd.DataFrame({"bases":(l_seq_res["pass_bases"] + l_seq_res["fail_bases"]).tolist(), "bases_line": l_seq_res["fail_bases"].tolist(), "classification": "passed", "colour": "#1F78B4", "left": l_seq_res["left"].tolist(), "right": l_seq_res["right"].tolist()})
            else:
                failed = pd.DataFrame({"reads":l_seq_res["fail_reads"].tolist(), "reads_line": 0, "classification": "failed", "colour": "#A6CEE3", "left": l_seq_res["left"].tolist(), "right": l_seq_res["right"].tolist()})
                passed = pd.DataFrame({"reads":(l_seq_res["pass_reads"] + l_seq_res["fail_reads"]).tolist(), "reads_line": l_seq_res["fail_reads"].tolist(), "classification": "passed", "colour": "#1F78B4", "left": l_seq_res["left"].tolist(), "right": l_seq_res["right"].tolist()})
            passed = passed.append(failed, sort=False)
        else:
            if normalised:
                passed = pd.DataFrame({"bases":l_seq_res["pass_bases"].tolist(), "bases_line": 0, "classification": "passed", "colour": "#1F78B4", "left": l_seq_res["left"].tolist(), "right": l_seq_res["right"].tolist()})
            else:
                passed = pd.DataFrame({"reads":l_seq_res["pass_reads"].tolist(), "reads_line": 0, "classification": "passed", "colour": "#1F78B4", "left": l_seq_res["left"].tolist(), "right": l_seq_res["right"].tolist()})
                
        if normalised:
            passed = passed.loc[passed.bases > 0,]
            plot_base = 'bases_line'
            plot_key = 'bases'
            plot_legend = "count (bases)"
        else:
            passed = passed.loc[passed.reads > 0,]
            logging.info("using read counts instead of bases!")
            plot_base = 'reads_line'
            plot_key = 'reads'
            plot_legend = "count (reads)"

        p = figure(title="Histogram showing read-length distribution", 
                   background_fill_color="lightgrey", plot_width=plot_width, 
                   plot_height=plot_height, tools=plot_tools,
                   x_range=Range1d(0, longest_read))
        
        p.quad(source=passed, top=plot_key, bottom=plot_base, left='left', right='right',
              fill_color='colour', line_color="white", legend_field='classification')
  
        if annotate_mean:
            vline = Span(location=mean_val, dimension='height', line_color='red', line_width=2)
            p.renderers.extend([vline])
            p.add_layout(LabelSet(x='x', y='y', text='text', level='glyph',
                                  source=ColumnDataSource(data=dict(x=[mean_val],
                                                                    y=[passed[plot_key].max()],
                                                                    text=['Mean'])),
                                  render_mode='canvas', text_align='right', text_color="red"))


        if annotate_n50:

            vline = Span(location=N, dimension='height', line_color='orange', line_width=2)
            p.renderers.extend([vline])
            p.add_layout(LabelSet(x='x', y='y', text='text', level='glyph',
                                  source=ColumnDataSource(data=dict(x=[N],
                                                                    y=[passed[plot_key].max()],
                                                                    text=['N50'])),
                                  render_mode='canvas', text_align='left', text_color="orange"))
        p.y_range.start = 0
        p.legend.location = "center_right"
        p.xaxis.axis_label = 'Sequence length (nt)'
        p.yaxis.axis_label = plot_legend
        p.yaxis.formatter = NumeralTickFormatter(format="0,0")
        p.xaxis.formatter = NumeralTickFormatter(format="0,0")
        p.grid.grid_line_color = "white"

        return self.handle_output(p, plot_type)   
    
    
    @functools.lru_cache()
    def extract_quality_stratified_data(self, bins=30, min=None, max=None):
    
        if min is None:
            min = self.seq_sum.mean_qscore_template.min()
        if max is None:
            max = self.seq_sum.mean_qscore_template.max()
        
        q_seq_sum = self.seq_sum[["mean_qscore_template", "passes_filtering"]]
        boundaries = np.linspace(min, max, 
                                 num=bins, endpoint=True, retstep=False)
        assignments = np.digitize(self.seq_sum.mean_qscore_template, boundaries)
        q_seq_sum = q_seq_sum.reindex(columns=q_seq_sum.columns.tolist()+["batch", "pass_reads", "fail_reads"])
        
        q_seq_sum.iloc[:,[2]] = assignments
        q_seq_sum.iloc[:,[3,4]] = 0
        q_seq_sum.loc[~q_seq_sum.passes_filtering, ["fail_reads"]] = 1
        q_seq_sum.loc[q_seq_sum.passes_filtering, ["pass_reads"]] = 1
        
        q_seq_res = q_seq_sum.groupby(["batch"]).agg(
            {"batch": "first",  "fail_reads": np.sum, "pass_reads": np.sum})

        # handle the possibility of missing values ... e.g. Dolores
        q_seq_res = q_seq_res.reindex(pd.Index(pd.Series(boundaries).index, name="hh")).reset_index()
        q_seq_res = q_seq_res.drop("hh", axis=1)
        # shape the data for quad presentation ...
        q_seq_res = q_seq_res.reindex(columns=q_seq_res.columns.tolist()+["left", "right"])
        q_seq_res.loc[1:,["left"]] = boundaries[:-1, np.newaxis]
        q_seq_res.loc[1:,["right"]] = boundaries[1:, np.newaxis]
        
        return boundaries, q_seq_res
    

    def plot_q_distribution(self, bins=30, **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        boundaries, q_seq_res = self.extract_quality_stratified_data(bins=bins)
        
        # and mung the data into something cleaner
        failed = pd.DataFrame({"reads":q_seq_res["fail_reads"], "reads_line": 0, "classification": "failed", "colour": "#A6CEE3", "left": q_seq_res["left"], "right": q_seq_res["right"]})
        passed = pd.DataFrame({"reads":(q_seq_res["pass_reads"] + q_seq_res["fail_reads"]), "reads_line": q_seq_res["fail_reads"], "classification": "passed", "colour": "#1F78B4", "left": q_seq_res["left"], "right": q_seq_res["right"]})
        passed = passed.append(failed, sort=False)
        
        plot_base = 'reads_line'
        plot_key = 'reads'
        plot_legend = "count (reads)"
        p = figure(title="Histogram showing distribution of quality values", 
                   background_fill_color="lightgrey", plot_width=plot_width,
                   plot_height=plot_height, tools=plot_tools)
        p.quad(source=passed, top=plot_key, bottom=plot_base, left='left', right='right',
               fill_color='colour', line_color="white", legend_field='classification', alpha=0.7)

        vline = Span(location=7, dimension='height', line_color='green', line_width=2)
        p.renderers.extend([vline])
        p.add_layout(LabelSet(x='x', y='y', text='text', level='glyph', 
                              source=ColumnDataSource(data=dict(x=[7], y=[passed['reads'].max()], text=['Q-filter'])),
                              render_mode='canvas', text_align='left', text_color="green"))

        p.y_range.start = 0
        p.legend.location = "center_right"
        p.xaxis.axis_label = 'Quality value (Phred)'
        p.yaxis.axis_label = plot_legend
        p.yaxis.formatter = NumeralTickFormatter(format="0,0")
        p.grid.grid_line_color = "white"

        return self.handle_output(p, plot_type)
    

    def plot_q_l_density(self, xbins=100, ybins=100, longest_read="Q0.99",
                         min_q="min", max_q="Q0.99", plot_depth_threshold=1, **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)

        if min_q is None:
            min_q = self.seq_sum.mean_qscore_template.min()
        if max_q is None:
            max_q = self.seq_sum.mean_qscore_template.max()
            
        min_q = self.get_length_character(
            key=min_q, passed=True, field="mean_qscore_template")
        max_q = self.get_length_character(
            key=max_q, passed=True, field="mean_qscore_template")
        
        longest_read = self.get_length_character(key="max", passed=True)
        shortest_read = self.get_length_character(key="min", passed=True)

        q_boundaries = np.linspace(min_q, max_q, num=ybins, endpoint=True, retstep=False)
        # l_boundaries = np.linspace(np.log10(100), np.log10(longest_read), num=xbins, endpoint=True)
        l_boundaries = np.logspace(np.log10(shortest_read), np.log10(longest_read), num=xbins)

        binned2d = stats.binned_statistic_2d(
            self.seq_sum['sequence_length_template'],
            self.seq_sum['mean_qscore_template'],
            np.repeat(1, len(self.seq_sum)), 'count',
            bins=[l_boundaries, q_boundaries]
        )
        layout = pd.DataFrame(binned2d.statistic)
        # layout = layout.transpose()
        layout.columns = binned2d.y_edge[:-1]
        layout.index = binned2d.x_edge[:-1]

        layout = layout.reset_index().melt(id_vars='index')
        layout.columns = ['column', 'row', 'count']
        layout = layout.fillna(0)

        # with the log scale we need a width column ...
        layout["width"] = layout["column"] / 10
        layout["alpha"] = 1
        layout.loc[layout["count"] <= plot_depth_threshold, 'alpha'] = 0

        logging.debug(layout)

        colors = []

        def colorFader(c1, c2, mix=0):  # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
            c1 = np.array(mpl.colors.to_rgb(c1))
            c2 = np.array(mpl.colors.to_rgb(c2))
            return mpl.colors.to_hex((1 - mix) * c1 + mix * c2)

        c1 = '#A6CEE3'
        c2 = 'blue'  
        n = 75

        for x in range(n + 1):
            colors.append(colorFader(c1, c2, x / n))

        mapper = LinearColorMapper(palette=colors, low=layout['count'].min() + 1, high=layout['count'].max())

        p = figure(title="Density plot showing relationship between quality and read length",
                   x_axis_location="below", plot_width=plot_width, plot_height=plot_height,
                   tools=plot_tools, toolbar_location='below', x_axis_type="log",
                   background_fill_color="lightgrey")

        p.rect(x="column", y="row", width="width", height=max_q / ybins,
               source=layout,
               fill_color={'field': 'count', 'transform': mapper},
               line_color=None,
               fill_alpha="alpha")

        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="10pt",
                             ticker=BasicTicker(desired_num_ticks=10),
                             # formatter=PrintfTickFormatter(format="%d%%"),
                             title="#reads",
                             label_standoff=6, border_line_color=None, location=(0, 0))
        p.add_layout(color_bar, 'right')

        hline = Span(location=7, dimension='width', line_color='green', line_width=2)
        p.renderers.extend([hline])
        p.add_layout(
            LabelSet(x='x', y='y', text='text', level='glyph', source=ColumnDataSource(data=dict(x=[longest_read],
                                                                                                 y=[7],
                                                                                                 text=['Q-filter'])),
                     render_mode='css', text_align='right', text_color="green"))

        vline = Span(location=self.get_length_character(key="mean", passed=True), dimension='height', line_color='red', line_width=2)
        p.renderers.extend([vline])
        p.add_layout(LabelSet(x='x', y='y', text='text', level='glyph', 
                              source=ColumnDataSource(data=dict(x=[self.get_length_character(key="mean", passed=True)], y=[max_q], text=['Mean'])),
                              render_mode='css', text_align='left', text_color="red"))
        p.xaxis.formatter = NumeralTickFormatter(format="0,0")
        p.xaxis.axis_label = 'Read length (nt)'
        p.yaxis.axis_label = 'Phred score (Q)'
        # show(p)
        return self.handle_output(p, plot_type)
    

    def plot_time_duty_reads(self, interval_mins=60, cumulative=True, 
                             include_total=False, include_failed=True, 
                             **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        (boundaries, temporal_seq_res) = self.extract_temporal_data(interval_mins)
        t_seq_res = temporal_seq_res.copy(deep=True)
        t_seq_res = t_seq_res.loc[t_seq_res.batch > 0,:]
        
        plot = figure(title='Plot showing sequence throughput against time', x_axis_label='Time (hours)',
                      y_axis_label='Sequence reads (n)', background_fill_color="lightgrey",
                      plot_width=plot_width, plot_height=plot_height, tools=plot_tools)
        plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
        
        if cumulative:
            if include_total:
                plot.line(boundaries[t_seq_res.index-1], np.cumsum(t_seq_res.counter), line_width=2, line_color='black', 
                          legend_label='Total reads')
            plot.line(boundaries[t_seq_res.index-1], np.cumsum(t_seq_res.pass_reads), line_width=2, line_color='#1F78B4',
                  legend_label='Passed reads')
            if include_failed:
                plot.line(boundaries[t_seq_res.index-1], np.cumsum(t_seq_res.fail_reads), line_width=2, line_color='#A6CEE3', 
                          legend_label='Failed reads')
            plot.legend.location = "top_left"
        else:
            if include_total:
                plot.line(boundaries[t_seq_res.index-1], t_seq_res.counter, line_width=2, line_color='black', 
                          legend_label='Total reads')
            plot.line(boundaries[t_seq_res.index-1], t_seq_res.pass_reads, line_width=2, line_color='#1F78B4',
                  legend_label='Passed reads')
            if include_failed:
                plot.line(boundaries[t_seq_res.index-1], t_seq_res.fail_reads, line_width=2, line_color='#A6CEE3', 
                          legend_label='Failed reads')
        return self.handle_output(plot, plot_type)


    @functools.lru_cache()
    def extract_temporal_data(self, interval_mins=60):
        boundaries = np.linspace(0, self.get_runtime(units='hours'),
                                 num=int(self.get_runtime(units='hours') * 60 / interval_mins + 1),
                                 endpoint=True, retstep=False)
        assignments = np.digitize(self.seq_sum['start_time'] / 60 / 60, boundaries)

        t_seq_sum = self.seq_sum[["sequence_length_template", "passes_filtering", "channel", "duration"]]
        
        t_seq_sum =  t_seq_sum.reindex(columns=t_seq_sum.columns.tolist() + 
                                       ["batch", "counter", "pass_bases",
                                        "fail_bases", "pass_reads", "fail_reads",
                                        "rate"])     
        t_seq_sum.loc[:,["batch"]] = assignments
        t_seq_sum.loc[:,["counter"]]=1
        t_seq_sum.loc[:,["pass_bases","fail_bases", "pass_reads", "fail_reads", "rate"]]=0
        
        t_seq_sum.loc[:, "rate"] = t_seq_sum['sequence_length_template'] / t_seq_sum['duration']
        
        t_seq_sum.loc[~t_seq_sum.passes_filtering, ["fail_reads"]] = 1
        t_seq_sum.loc[t_seq_sum.passes_filtering, ["pass_reads"]] = 1
        t_seq_sum.loc[t_seq_sum.passes_filtering, ["pass_bases"]] = t_seq_sum.loc[t_seq_sum.passes_filtering, ["sequence_length_template"]].iloc[:,0]
        t_seq_sum.loc[~t_seq_sum.passes_filtering, ["fail_bases"]] = t_seq_sum.loc[~t_seq_sum.passes_filtering, ["sequence_length_template"]].iloc[:,0]
        
        t_seq_res = t_seq_sum.groupby(["batch"]).agg(
            {"batch": "first", "sequence_length_template": np.sum, 
             "counter": np.sum, "fail_reads": np.sum, "pass_reads": np.sum,
             "pass_bases": np.sum, "fail_bases": np.sum, 
             "channel": pd.Series.nunique, "rate": 
                 lambda x: (x.quantile(q=0.25), x.quantile(q=0.5), x.quantile(q=0.75))})
        t_seq_res = t_seq_res.reindex(pd.Index(pd.Series(boundaries).index, name="hh")).reset_index()
        t_seq_res = t_seq_res.fillna(0)
        # remove 0 val - this is added during feature filling ...
        t_seq_res = t_seq_res.loc[t_seq_res.index > 0,:]
        return (boundaries, t_seq_res.drop("hh", axis=1))
        

    def plot_time_duty_bases(self, interval_mins=60, scale="Gigabases", 
                             cumulative=True, milestones=[0.5, 0.9], 
                             include_total=False, include_failed=True,
                             **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        # validate the scale ...
        scaleVal = 1
        if scale == "Gigabases":
            scaleVal = 1e9
        elif scale == "Megabases":
            scaleVal = 1e6
        elif scale == "Kilobases":
            scaleVal = 1e3

        (boundaries, temporal_seq_res) = self.extract_temporal_data(interval_mins)
        t_seq_res = temporal_seq_res.copy(deep=True)
        t_seq_res = t_seq_res.loc[t_seq_res.batch > 0,:]
        
        # this deep data copy is required because we are scaling the data
        # and this messes with the core cached object ...
        t_seq_res.iloc[:,[1,5,6]] = t_seq_res.iloc[:,[1,5,6]] / scaleVal
        
        plot = figure(title='Plot showing sequence throughput against time', x_axis_label='Time (hours)',
                      y_axis_label='Sequence {} (n)'.format(scale), background_fill_color="lightgrey",
                      plot_width=plot_width, plot_height=plot_height, tools=plot_tools)
        if cumulative:
            if include_total:
                plot.line(boundaries[t_seq_res.index-1], np.cumsum(t_seq_res.sequence_length_template), line_width=2, line_color='black', 
                          legend_label='bases across all reads')
            plot.line(boundaries[t_seq_res.index-1], np.cumsum(t_seq_res.pass_bases), line_width=2, line_color='#1F78B4',
                  legend_label='bases from passed reads')
            if include_failed:
                plot.line(boundaries[t_seq_res.index-1], np.cumsum(t_seq_res.fail_bases), line_width=2, line_color='#A6CEE3', 
                          legend_label='bases from failed reads')
            for milestone in milestones:
                (bases, times) = self.get_sequence_base_point(fraction=milestone, interval_mins=interval_mins)
                bases = bases / scaleVal
                legend = "T{:.0f}".format(milestone * 100)
                plot.line([0, times, times], [bases, bases, 0], line_width=2, line_color='red')
                plot.add_layout(Label(x=times, y=bases, text=legend, text_color='red'))
            plot.legend.location = "top_left"
        else:
            if include_total:
                plot.line(boundaries[t_seq_res.index-1], t_seq_res.sequence_length_template, line_width=2, line_color='black', 
                          legend_label='bases across all reads')
            plot.line(boundaries[t_seq_res.index-1], t_seq_res.pass_bases, line_width=2, line_color='#1F78B4',
                  legend_label='bases from passed reads')
            if include_failed:
                plot.line(boundaries[t_seq_res.index-1], t_seq_res.fail_bases, line_width=2, line_color='#A6CEE3', 
                          legend_label='bases from failed reads')
        return self.handle_output(plot, plot_type)


    def get_sequence_base_point(self, fraction=0.5, interval_mins=5, column="pass_bases"):
        (boundaries, t_seq_res) = self.extract_temporal_data(interval_mins)
        target_value = t_seq_res[column].sum() * fraction
        return (target_value, np.interp(target_value, np.cumsum(t_seq_res[column]), boundaries[:-1]))



    def plot_translocation_speed(self, interval_mins=60, **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        
        (boundaries, temporal_seq_res) = self.extract_temporal_data(interval_mins)
        t_seq_res = temporal_seq_res.copy(deep=True)
        t_seq_res = t_seq_res.loc[t_seq_res.batch > 0]

        iqrd = pd.DataFrame(t_seq_res.rate.tolist(), columns=["q1","q2","q3"])
        iqrd["iqr"] = iqrd["q3"] - iqrd["q1"]
        
        iqrd["upper"] = iqrd["q3"] + 1.5 * iqrd["iqr"]
        iqrd["lower"] = iqrd["q1"] - 1.5 * iqrd["iqr"]

        plot = figure(title='Plot showing sequencing rate against time', x_axis_label='Time (hours)',
                      y_axis_label='Sequencing rate (bases/s)', background_fill_color="lightgrey",
                      plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        time_points = boundaries[t_seq_res.index.tolist()]

        plot.segment(time_points, iqrd["upper"], time_points, iqrd["q3"], line_color="black")
        plot.segment(time_points, iqrd["lower"], time_points, iqrd["q1"], line_color="black")
        plot.vbar(time_points, 0.7, iqrd["q2"], iqrd["q3"], fill_color="#E08E79", line_color="black")
        plot.vbar(time_points, 0.7, iqrd["q1"], iqrd["q2"], fill_color="#3B8686", line_color="black")

        plot.rect(time_points, iqrd["lower"], 0.2, 0.01, line_color="black")
        plot.rect(time_points, iqrd["upper"], 0.2, 0.01, line_color="black")

        return self.handle_output(plot, plot_type)

    def plot_functional_channels(self, interval_mins=60, **kwargs):
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        (boundaries, temporal_seq_res) = self.extract_temporal_data(interval_mins)
        t_seq_res = temporal_seq_res.copy(deep=True)
        t_seq_res = t_seq_res.loc[t_seq_res.batch > 0]

        time_points = boundaries[t_seq_res.index.tolist()]

        plot = figure(title='Plot showing number of observed channels against time', x_axis_label='Time (hours)',
                      y_axis_label='Number of active channels (n)', background_fill_color="lightgrey",
                      plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        plot.step(time_points, t_seq_res.channel, line_width=2, mode="before")
        return self.handle_output(plot, plot_type)

    def is_barcoded_dataset(self):
        """
        There are a couple of different ways that we need to consider for
        barcode context - has the dataset been demiultiplexed at time of
        basecalling within e.g. Guppy - or as a downstream process
        
        key information relied upon here is the barcode field in the 
        sequencing summary file ...

        Returns
        -------
        bool
            DESCRIPTION.

        """
        return "barcode_arrangement" in self.seq_sum.columns


    def plot_barcode_info(self, threshold=1, **kwargs):
        if not self.is_barcoded_dataset():
             return None
        (plot_width, plot_dpi) = self.handle_kwargs(["plot_width", "plot_dpi"], **kwargs)

        barcodes = self.tabulate_barcodes(threshold=threshold)
        bi = pd.Series(barcodes.index)

        classified = barcodes.iloc[bi[bi != "unclassified"].index.values]
        unclassified = barcodes.iloc[bi[bi == "unclassified"].index.values]

        bc_fraction = InfographicNode(
            legend="Reads with barcode",
            value="{:.2f} %".format(
                classified['count'].sum() / 
                (classified['count'].sum() + 
                 unclassified['count'].sum()) * 100),
            graphic='chart-pie')
        bc_bc_count = InfographicNode(
            legend="Barcoded libraries",
            value="{}".format(len(classified.index)), 
            graphic='barcode')
        bc_variance = InfographicNode(
            legend="Barcode variance",
            value=">{}\n{}<".format(
                classified["count"].min(),classified["count"].max()),
            graphic='sort-numeric-down')
        infographic_data = [bc_fraction, bc_bc_count, bc_variance]
        ip = InfographicPlot(infographic_data, rows=1, columns=3)
        return ip.plot_infographic(plot_width, plot_dpi)


    def calculate_mean_quality(self, series):
        series = series.dropna()
        return -10 * log10((10 ** (series / -10)).mean())

    @functools.lru_cache()
    def tabulate_barcodes(self, threshold=1):
        if not self.is_barcoded_dataset():
             return None
        
        bc_seq_sum = self.seq_sum.loc[
            self.seq_sum.passes_filtering, [
                "sequence_length_template", "barcode_arrangement", 
                "mean_qscore_template"]]
        bc_seq_sum["count"]=1
        total_reads = np.sum(bc_seq_sum["count"])
        
        def get_perc(x):
            return x.sum()/total_reads * 100
            
        bc_res = bc_seq_sum.groupby(
            "barcode_arrangement"
            ).agg({"count": [np.sum, get_perc],
                   "mean_qscore_template": [self.calculate_mean_quality], 
                   "sequence_length_template":[np.sum, np.min, np.max, 
                                               np.mean, self.calculate_n_val]})
        bc_res.columns = bc_res.columns.droplevel()
        bc_res.index.name = ''
        bc_res.columns = ["count", "%", "mean_q", "Mbases", "min", "max", 
                          "mean", "N50"]
        
        bc_res["mean_q"] = bc_res["mean_q"].round(2)
        bc_res["Mbases"] = bc_res["Mbases"].round(2)
        bc_res["mean"] = bc_res["mean"].round(2)
        bc_res["%"] = bc_res["%"].round(2)
        
        bc_res = bc_res.loc[bc_res["%"] > threshold,:]
        
        return bc_res

    def plot_barcodes(self, threshold=1, **kwargs):
        (
            plot_width, plot_height, plot_type, plot_tools
         ) = self.handle_kwargs(
             ["plot_width", "plot_height", "plot_type", "plot_tools"], 
             **kwargs)
        
        if not self.is_barcoded_dataset():
            return None
        barcodes = self.tabulate_barcodes(threshold=threshold)
         
        p = figure(title="Histogram showing abundance of different barcodes", 
                   background_fill_color="lightgrey", plot_width=plot_width,
                   plot_height=plot_height, tools=plot_tools,
                   x_range=barcodes.index.tolist())

        p.vbar(barcodes.index.tolist(), top=barcodes["count"].tolist(),
               width=0.75, fill_alpha=0.7, fill_color="#1F78B4")

        p.y_range.start = 0
        p.xaxis.axis_label = 'Barcode assignment'
        p.yaxis.axis_label = "Number of Reads"
        p.yaxis.formatter = NumeralTickFormatter(format="0,0")
        p.xaxis.major_label_orientation = math.pi/2
        p.grid.grid_line_color = "white"
        
        return self.handle_output(p, plot_type)
    
    
    def barcode_subset(self, barcode_id):
        if not self.is_barcoded_dataset():
             return None
        barcodes = self.tabulate_barcodes()
        if barcode_id not in barcodes.index:
            logging.warning("requested barcode not present in barcode file")
            return None
        
        subset = self.seq_sum[self.seq_sum["barcode_arrangement"]==barcode_id]
        sub_ssh = SequenceSummaryHandler(
            target_data=subset, fcid="{}\n{}".format(
                self.get_flowcell_id(), barcode_id))
        self.sync(sub_ssh)
        return sub_ssh





class SequencingSummaryGetChannelMap:
    """
    This class is responsible for the handling of flowcell channel maps -
    the prototype of this class was written in R (see @sagrudd/nanopoRe);
    code has been transposed and simplified
    """

    def __init__(self, seq_sum):
        self.seq_sum = seq_sum

    def get_platform(self):
        """
        method scores the defined channels in the provided sequencing_summary
        to look for the largest defined channel - based on the number observed
        reports whether it is most likely to be Flongle / MinION or PromethION

        Returns
        -------
        String representation of flowcell type (MinION/Flongle/PromethION)

        """

        platform = "MinION"
        max_channel = self.seq_sum['channel'].max()
        logging.debug("MaxChannel == ", max_channel)

        if max_channel < 130:
            # this is likely to be a Flongle ...
            platform = "Flongle"

        if max_channel > 1000:
            # this is likely to be a PromethION
            platform = "PromethION"
        logging.debug("flowcell_type identified as %s" % platform)
        return platform

    def get_minion_map(self):
        """
        The R code is below; straight forward and minimal ...
        
        https://wiki/pages/viewpage.action?spaceKey=ELEC&title=Minion+Chip+map
        
        blockCalc <- function(i) {
            m <- matrix(seq(i, i + 63, by = 1), ncol = 8, byrow = TRUE)
            cbind(m[seq(5, 8, by = 1), ], m[seq(4), rev(seq(8))])
        }
        layout <- do.call(rbind, lapply(
            c(1, 449, 385, 321, 257, 193, 129, 65), blockCalc))
        # transpose the layout for cleaner presentation ...
        layout <- t(layout)

        channelMap <- as.data.frame(cbind(channel=as.vector(t(layout)),which(
        layout == as.vector(layout), arr.ind = TRUE)))

        Returns
        -------
        None.

        """

        def block_calc(i):
            m = np.arange(i, (i + 64)).reshape((8, 8), order='C')
            row = np.c_[m[np.arange(4, 8).tolist(), :],
                        m[np.arange(0, 4).tolist(),][:, np.arange(8)[::-1]]
            ]
            return row

        vector = [1, 449, 385, 321, 257, 193, 129, 65]
        coord_vector = pd.Series(vector).apply(block_calc)
        layout = np.vstack(np.stack(coord_vector))
        coords_df = pd.DataFrame(layout).reset_index().melt(id_vars='index')
        coords_df.columns = ['row', 'column', 'channel']
        return coords_df

    def get_flongle_map(self):
        """
        This again derived from the nanopoRe package implemented previously
        R code contained here for reference

        layout <- matrix(c(seq(1, 12), 0, seq(13, 24), 0, seq(25, 114), 0,
                           seq(115, 126), 0), ncol = 13, byrow = TRUE)
        layout <- layout[rev(seq(10)), ]
        channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)),
            which(layout == as.vector(layout), arr.ind = TRUE)))

        Returns
        -------
        None.

        """
        layout = np.concatenate([np.arange(1, 13),
                                 np.array([0]),
                                 np.arange(13, 25),
                                 np.array([0]),
                                 np.arange(25, 115),
                                 np.array([0]),
                                 np.arange(115, 127),
                                 np.array([0])]).reshape(10, 13)
        coords_df = pd.DataFrame(layout).reset_index().melt(id_vars='index')
        coords_df.columns = ['row', 'column', 'channel']
        return coords_df

    def get_promethion_map(self):
        """
        
        chunk <- function(i) {
            m <- matrix(seq_len(250), ncol=10, byrow=TRUE)
            m + i
        }
        layout <- do.call(cbind, lapply(seq(from=0, to=2750, by=250), chunk))
        channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)),
        which(layout == as.vector(layout), arr.ind = TRUE)))

        Returns
        -------
        None.

        """

        def chunk(i):
            return np.arange(1, 251).reshape(25, 10) + i

        layout = np.hstack(
            np.stack(pd.Series(np.arange(0, 2751, 250)).apply(chunk)))
        coords_df = pd.DataFrame(layout).reset_index().melt(id_vars='index')
        coords_df.columns = ['row', 'column', 'channel']
        return coords_df

    def get_platform_map(self):
        """
        simplifies the selection of the most appropriate channel layout map

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        platform = self.get_platform()
        defined_platforms = {"MinION": self.get_minion_map(),
                             "Flongle": self.get_flongle_map(),
                             "PromethION": self.get_promethion_map()
                             }
        return defined_platforms[platform]

    def get_platform_density(self):
        platform_map = self.get_platform_map()
        # collapse the channel data to get count information ...
        channel_counts = self.seq_sum.groupby(
            'channel')['channel'].agg(['count'])
        # merge the count data with the map coordinates
        channel_coords_counts = pd.merge(
            platform_map, channel_counts, on='channel', how='left')
        if self.get_platform() == "MinION":
            logging.info("rotating MinION data")
            channel_coords_counts = channel_coords_counts.rename(
                columns={"row": "column", "column": "row"})
        return channel_coords_counts
    
    
class FixSequencingSummary:
    
    def __init__(self, fin):
        self.fin = fin
        
    def parse(self):
        lines = 0
        fh = open(self.fin, 'r')
        while 1:
            line = fh.readline()
            if not line:
                break
            line = line.strip()
            fields = len(line.split("\t"))
            print(fields)
            lines += 1
            if lines > 10: 
                fh.close()
                return fields
                
    def parse_fields(self, count):
        fh = open(self.fin, 'r')
        print("filesize ?? %s" % os.fstat(fh.fileno()).st_size)

        while 1:
       
            line = fh.readline()
            if not line:
                break
            line = line.strip()
            fields = len(line.split("\t"))
            if fields != count:
                print("%s fields at %s" % (fields, fh.tell()))
 
        fh.close()
                
    def fixit(self, count):
        fh = open(self.fin, 'r')
        fw = open(self.fin+".fixed", 'w')
        with tqdm(ascii=True, total=os.fstat(fh.fileno()).st_size) as pbar:
            while 1:
                pbar.update(fh.tell())
                line = fh.readline()
                if not line:
                    break
                line = line.strip()
                fields = len(line.split("\t"))
                if fields == count:
                    fw.write(line+"\n")
        pbar.close()
        fw.close()
                
        
                
                
