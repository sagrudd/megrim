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
from genome_geometry import GenomeGeometry
from bisect import bisect_left
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from time import time
from infographic_plots import InfographicPlot, InfographicNode
from bokeh.io import export_png, show
from bokeh.models import LinearColorMapper, BasicTicker, PrintfTickFormatter, ColorBar, Label, LabelSet, NumeralTickFormatter
from bokeh.plotting import figure
from palettable.colorbrewer.sequential import Blues_9


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
        


class SequenceSummaryHandler:
    
    def __init__(self, target_file):
        self.target_file = target_file
        self._import_data()
        
    def _import_data(self):
        """
        method parses the provided Sequencing_summary.txt file and loads 
        the reduced contents into memory - a DASK dataframe is used to
        allow for the scaling to large PromethION runs (or equivalent)

        Returns
        -------
        None.

        """
        self.seq_sum = dd.read_csv(
            self.target_file, 
            delimiter='\t',
            blocksize=64000000,
        )
        
        # slice out the head entries for further functionality
        self.seq_sum_head = self.seq_sum.head()
        
        # start excluding dask columns that are not of core interest
        keep = ['channel', 'start_time', 'duration', 'num_events', 
                'sequence_length_template', 'mean_qscore_template', 
                'passes_filtering',
        ]
        for col in self.seq_sum.columns:
            if not col in keep:
                print("dropping %s" % col)
                self.seq_sum = self.seq_sum.drop(col, axis=1)
        pbar = ProgressBar()
        pbar.register()
        with Timer("Elapsed time to extract sequence data {}"):
                self.seq_sum = self.seq_sum.persist()
        pbar.unregister()
        
        
    def get_flowcell_id(self):
        """
        This method pulls the flowcell id from the fastq filename of the 
        original sequence ... this code may need to be maintained ...

        Returns
        -------
        fastq_id : TYPE
            DESCRIPTION.

        """
        fastq_id = str(self.seq_sum_head['filename_fastq'][0]).split("_")
        if len(fastq_id) == 3:
            fastq_id = fastq_id[0]
        logging.debug("flowcell read as [%s]" % fastq_id)
        #print(("_", fastq_id))
        return fastq_id
        

    def executive_summary(self):
        """
        Method prepares a three panel infographic plot that summarises the
        flowcell, reads and bases sequenced within the associated run

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        read_count = len(self.seq_sum)
        total_bases = self.seq_sum['sequence_length_template'].sum().compute()
        # passed_bases = self.seq_sum[self.seq_sum['passes_filtering']]['sequence_length_template'].sum().compute()
        
        
        flowcell_node = InfographicNode(legend="flowcell", 
                                value=self.get_flowcell_id(), 
                                graphic='fingerprint')

        readcount_node = InfographicNode(legend="Reads produced", 
                                value="{:,}".format(read_count),
                                graphic='filter')

        gb_seq_val = InfographicNode(legend="Gigabases called", 
                                value="{:.2f}".format(total_bases/1e9),
                                graphic='flag-checkered')
        infographic_data = [flowcell_node, readcount_node, gb_seq_val]
        ip = InfographicPlot(infographic_data, rows=1, columns=3)  
        return ip.plot_infographic() 

    def get_absolute_runtime(self):
        max_time = self.seq_sum['start_time'].max().compute()
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
                 "hours": 60*60,
                 "days": 60*60*24}
        canonical_runs = [1,2,4,6,8,12,18,24,36,48,60,72]
        
        if not units in scale.keys():
            raise Exception("{%s} is not a qualified runtime unit ...")
        
        runtime = math.ceil(self.get_absolute_runtime() / scale[units])        
        if rounded:
            def take_closest(values, val):
                pos=bisect_left(values, val)
                if pos == len(values):
                    return values[-1]
                else:
                    return values[pos]
            runtime = take_closest(canonical_runs, runtime)            
        return runtime
    
    
    
    def plot_passed_gauge(self):
        read_count = len(self.seq_sum)
        passed_read_count = self.seq_sum.passes_filtering.sum().compute()
        perc_val = passed_read_count / read_count * 100
        
        p = figure(plot_width=600, plot_height=400, x_range=(0.25, 1.75), y_range=(0.7, 1.5), toolbar_location=None)
        
        start_val = 0
        middle_val = (math.pi / 100) * (100 - perc_val)
        end_val = math.pi
        
        p.annular_wedge(x=[1], y=[1], inner_radius=0.2, outer_radius=0.5,
                        start_angle=middle_val, end_angle=end_val, color="green", alpha=0.6)
        
        p.annular_wedge(x=[1], y=[1], inner_radius=0.2, outer_radius=0.5,
                        start_angle=start_val, end_angle=middle_val, color="orange", alpha=0.6)
        
        label = Label(x=1, y=1, text="{:.1f}%".format(perc_val), x_units='data', y_units='data', 
                      text_align='center', text_font_style='bold', text_font_size='1.5em')
        legend = Label(x=1, y=0.9, 
                       text="Percentage of reads passing QC filter", 
                       x_units='data', y_units='data', 
                      text_align='center', text_font_size='1.9em')
        p.add_layout(label)
        p.add_layout(legend)
        p.axis.visible = False
        p.xgrid.visible=False
        p.ygrid.visible=False
        show(p)
        
        return "ThisIsAFilename.png"

    def plot_channel_activity(self):
        channel_map = SequencingSummaryGetChannelMap(self.seq_sum)
        #layout = channel_map.get_platform_map()
        layout = channel_map.get_platform_density()
        layout = layout.fillna(0)
        layout['row'] = layout['row'].astype(str)
        layout['column'] = layout['column'].astype(str)
        layout['count'] = layout['count'].astype(int)
        
        print(layout)

        #colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
        colors = Blues_9.hex_colors
        mapper = LinearColorMapper(palette=colors, low=layout['count'].min(), high=layout['count'].max())

        TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

        rows = list(layout.row.unique())
        columns = list(layout.column.unique())

        p = figure(title="channel activity plot",
            x_range=columns, y_range=rows,
            x_axis_location="above", plot_width=1200, plot_height=800,
            tools=TOOLS, toolbar_location='below')

        p.axis.visible = False
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "5pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = math.pi / 3
        p.title.text_font_size = '18pt'

        p.rect(x="column", y="row", width=1, height=1,
               source=layout,
               fill_color={'field': 'count', 'transform': mapper},
               line_color=None)

        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="10pt",
                     ticker=BasicTicker(desired_num_ticks=len(colors)),
                     #formatter=PrintfTickFormatter(format="%d%%"),
                     title="#reads",
                     label_standoff=6, border_line_color=None, location=(0, 0))
        p.add_layout(color_bar, 'right')
        #show(p)
        export_png(p, filename="plot.png")
        return "ThisIsAFilename.png"
    
    def library_characteristics_infographic(self):
        
        geometry = GenomeGeometry(pd.Series(self.seq_sum[self.seq_sum['passes_filtering']]['sequence_length_template'].compute()))
        longest_read = geometry.get_longest_read()
        logging.debug("longest_read == %s" % (longest_read))
        mean_read_length = geometry.get_mean_length()
        logging.debug("mean_read_length == %s" % (mean_read_length))
        read_n50_length = geometry.get_n_value(n=50)
        logging.debug("read_n50_length == %s" % (read_n50_length))
        passed_mean_q = geometry.calculate_mean_quality(pd.Series(self.seq_sum[self.seq_sum['passes_filtering']]['mean_qscore_template'].compute()))
        logging.debug("passed_mean_q == %s" % (passed_mean_q))
        failed_mean_q = geometry.calculate_mean_quality(pd.Series(self.seq_sum[~self.seq_sum['passes_filtering']]['mean_qscore_template'].compute()))
        logging.debug("failed_mean_q == %s" % (failed_mean_q))
        
        #df$info <- c(passedMeanLength, N50, passedMeanQ, failedMeanQ, prettyNum(max(passedSeqs$sequence_length_template), big.mark=","))
        #df$key <- c("Mean Read Length (nt)","N50","Mean Read Quality (QV)","Mean Failed QV","Longest Read")
        #df$icon <- fontawesome(c("fa-bar-chart", "fa-play", "fa-area-chart", "fa-bug", "fa-sort"))
        
        mean_read_length_node = InfographicNode(legend="Mean Read Length (nt)", 
                                value="{:.2f}".format(mean_read_length), 
                                graphic='map-signs')
        read_n50_length_node = InfographicNode(legend="N50", 
                                value="{:,}".format(read_n50_length), 
                                graphic='bullseye')
        passed_mean_q_node = InfographicNode(legend="Mean Read Quality (QV)", 
                                value="{:.2f}".format(passed_mean_q), 
                                graphic='award')
        failed_mean_q_node = InfographicNode(legend="Mean Failed QV", 
                                value="{:.2f}".format(failed_mean_q), 
                                graphic='bug')        
        longest_read_node = InfographicNode(legend="Longest Read", 
                                value="{:,}".format(longest_read), 
                                graphic='sort')
        infographic_data = [mean_read_length_node, read_n50_length_node,
                            passed_mean_q_node, failed_mean_q_node,
                            longest_read_node]
        ip = InfographicPlot(infographic_data, rows=1, columns=5)  
        return ip.plot_infographic() 
    
    
    def plot_sequence_length(self, normalised=True, include_failed=True, bins=30):        
        # there are some approaches such as np.histogram; seems to split
        # data into clear bins; but not sure on how for stacked ranges ...
        # let's use a few more lines of code and perform manually
        
        geometry = GenomeGeometry(pd.Series(self.seq_sum[self.seq_sum['passes_filtering']]['sequence_length_template'].compute()))
        geometryF = GenomeGeometry(pd.Series(self.seq_sum[~self.seq_sum['passes_filtering']]['sequence_length_template'].compute()))
        
        longest_read = geometry.get_longest_read()
        boundaries = np.linspace(0, longest_read, num=bins, endpoint=True, retstep=False)
        print(boundaries)
        indsP = np.digitize(geometry.get_lengths(), boundaries)
        indsF = np.digitize(geometryF.get_lengths(), boundaries)
        countsP = np.unique(indsP, return_counts=True, return_inverse=True)
        countsF = np.unique(indsF, return_counts=True, return_inverse=True)
        print(countsP[1])
        
        
        def count_bases(x, assignments, reads):
            return reads[assignments==x].sum()
        
        chunksP = pd.Series(np.unique(countsP[1]))
        basesP = chunksP.apply(count_bases, assignments=countsP[1], reads=geometry.get_lengths())

        chunksF = pd.Series(np.unique(countsF[1]))
        basesF = chunksF.apply(count_bases, assignments=countsF[1], reads=geometryF.get_lengths())
        
        # pad the missing values from these ranges 
        dfP = pd.DataFrame({'bases_line': 0,
                            'count_line': 0,
                            'count': np.repeat(0, bins-1), 
                            'bases': np.repeat(0, bins-1), 
                            'classification': 'passed',
                            'colour': "#1F78B4",
                            'left': boundaries[:-1],
                            'right': boundaries[1:]})
        dfP.loc[bins-1] = np.array([0, 0, 0, 0, 'passed',"#1F78B4",dfP.right[bins-2],dfP.right[bins-2]+(dfP.right[bins-2]-dfP.left[bins-2])])
        dfP.loc[countsP[0]-1, 'count'] = countsP[2] 
        dfP.loc[countsP[0]-1, 'bases'] = basesP
        
        if include_failed:
            dfF = pd.DataFrame({'bases_line': 0,
                                'count_line': 0,
                                'count': np.repeat(0, bins-1), 
                                'bases': np.repeat(0, bins-1), 
                                'classification': 'failed',
                                'colour': "#A6CEE3",
                                'left': boundaries[:-1],
                                'right': boundaries[1:]})
            dfF.loc[bins-1]  = np.array([0, 0, 0, 0, 'failed',"#A6CEE3",dfF.right[bins-2],dfF.right[bins-2]+(dfF.right[bins-2]-dfF.left[bins-2])])
            dfF.loc[countsF[0]-1, 'count'] = countsF[2] 
            dfF.loc[countsF[0]-1, 'bases'] = basesF
            
            # one challenge here is that bokeh quads do not support stacking ...
            # this should be managed by self ...
            dfP['bases_line'] = dfF['bases']
            dfP['bases'] = dfP['bases_line']+ dfP['bases']
            dfP['count_line'] = dfF['count']
            dfP['count'] = dfP['count_line']+ dfP['count']
            dfP = dfP.append(dfF)
        
        print(dfP)
        
        plot_base = 'bases_line'
        plot_key = 'bases'
        plot_legend = "count (bases)"
        if not normalised:
            print("using counts instead of bases!")
            plot_base = 'count_line'
            plot_key = 'count'
            plot_legend = "count (reads)"
        
        p = figure(title="Histogram showing read-length distribution", background_fill_color="lightgrey")
        p.quad(source=dfP, top=plot_key, bottom=plot_base, left='left', right='right',
           fill_color='colour', line_color="white", legend_field= 'classification')
        
        
        p.y_range.start = 0
        p.legend.location = "center_right"
        p.xaxis.axis_label = 'Sequence length (nt)'
        p.yaxis.axis_label = plot_legend
        p.yaxis.formatter=NumeralTickFormatter(format="0,0")
        p.grid.grid_line_color="white"
        
        export_png(p, filename="plot3.png")
        return "ThisIsAFilename.png"
        
    

    
    
    
    def plot_q_distribution(self):
        i = 1
        return "ThisIsAFilename.png"
        
    def plot_q_l_density(self):
        i = 1
        return "ThisIsAFilename.png"
        
    def plot_time_duty(self, interval=0.25):
        i = 1
        return "ThisIsAFilename.png"
        
    def plot_cumulative_bases(self, interval=0.25, milestones=[0.5, 0.9]):
        i = 1
        return "ThisIsAFilename.png"
        
    def plot_cumulative_reads(self, interval=0.25, milestones=[0.5]):
        i = 1
        return "ThisIsAFilename.png"
    
    def plot_translocation_speed(self, interval=1):
        i = 1
        return "ThisIsAFilename.png"
    
    def plot_functional_channels(self, interval=1):
        i = 1
        return "ThisIsAFilename.png"
    
    def is_barcoded_dataset(self):
        return False
    
    def merge_barcode_summary(self, barcode_summary_file):
        i = 1
        
    def plot_barcode_info(self):
        i = 1
        return "ThisIsAFilename.png"
    
    def tabulate_barcodes(self, threshold=0.01):
        i = 1
        
    def plot_barcodes(self, threshold=0.01):
        data = self.tabulate_barcodes(threshold=threshold)
        return "ThisIsAFilename.png"
        

nul = """
ssh = SequenceSummaryHandler("/Users/srudd/Desktop/Dolores_PromethION_1M.txt")
ssh.executive_summary()
"""


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
            row = np.c_[m[np.arange(4, 8).tolist(),:],
                        m[np.arange(0, 4).tolist(),][:,np.arange(8)[: :-1]]
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
                                 np.array([0])]).reshape(10,13)
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
            return np.arange(1,251).reshape(25,10) + i
        layout = np.hstack(np.stack(pd.Series(np.arange(0, 2751, 250)).apply(chunk)))
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
        channel_counts = self.seq_sum.groupby('channel')['channel'].agg(['count']).compute()
        # merge the count data with the map coordinates
        return pd.merge(platform_map, channel_counts, on='channel', how='left')
        
        
    

        

        
    