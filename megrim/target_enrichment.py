#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
target_enrichment.py
====================

This module contains a single class for managing the core of the Cas9-based
target enrichment tutorial. The essence of the workflow is the judicious
application of pyranges and pysam to distil coverage information from
across the genome - emphasis on target regions defined within a starting
bed file.

please see the Nanopore Cas9 tutorial at
http://github.com/nanoporetech/ont_tutorial_cas9

Created on Tue Jan 21 18:09:10 2020
@author: srudd
"""


from megrim.environment import Flounder
from megrim.genome_geometry import BamHandler, BedHandler
from megrim.reference_genome import ReferenceGenome, augment_annotation, \
    weighted_percentile
from megrim.infographic_plots import InfographicPlot, InfographicNode
import numpy as np
import pandas as pd
import logging
import functools
from bokeh.plotting import figure
from bokeh.models import Span, NumeralTickFormatter


class TargetEnrichment(Flounder):
    r"""
    TargetEnrichment class for supporting the Cas9 tutorial.

    This class provides the target_enrichment functionality used within the
    ont_tutorial_cas9 tutorial. The class leverages functionality from other
    modules within the megrim framework to simplify and abstract the process
    of defining on-target regions of the genome, defining the background
    genome and performing calculations to derive off-target regions.

    Parameters
    ----------
        reference_genome: ``file.path``
            a reference genome is required for the analyses - the provided
            path will be used to create an instance of
            :class:`~megrim.reference_genome`
        bed_file: ``**file.path**``
            a bed format file is required to define the genomic locations
            of the target region(s) of interest
        bam_file: ``**file.path**``
            The bam file is required to provide mapping information for reads
            across the genome and within the target regions. Please note that
            this bam file *must* be indexed
        target_proximity: ``int``
            This defines the regions up- and down-stream of the defined
            target regions that will be considered as target-proximal
            regions. The default is 2000 bases
        background_threshold: ``**int**``
            This threshold is used to defined what is off-target - the
            mean depth-of-coverage across the non-target/target-proximal
            regions of the genome is calculated and off-target is defined
            as a region with > (background_threshold * mean background
            coverage). The default value is 20.
        tile_size: ``**int**``
            The genome is tiled for the depth of coverage analyses. This
            parameter is used to define the width of the tile. Shorter windows
            provide higher resolution but have a compute penalty. The default
            value is 1000.

    .. note::

        the TargetEnrichment class inherits the Flounder class and \
        requires functionality from classes including Bam, Bed and \
        ReferenceGenome
    """

    def __init__(self, reference_genome=None, bed_file=None, bam_file=None,
                 target_proximity=2000, background_threshold=20,
                 tile_size=1000):
        Flounder.__init__(self)

        logging.info("Preparing reference genome")
        self.ref = ReferenceGenome(reference_genome)

        logging.info("Preparing BED coordinates")
        self.bed = BedHandler(bed_file)
        self.bed.set_reference(self.ref)
        self.bed.set_target_proximity(target_proximity)

        logging.info("setting BAM file")
        self.bam = BamHandler(bam_file)
        self.background_threshold = background_threshold
        self.target_proximity = target_proximity
        self.tile_size = tile_size

    @functools.lru_cache()
    def get_on_target(self):
        """
        Get pyranges format on-target coordinates and coverage.

        This method returns the coordinates for the defined target regions
        used within the target enrichment analysis. This is a simple
        reflection of the content provided within the provided BED file

        Returns
        -------
        on_target_universe : pyranges
            This is a pyranges refection of the tab-delimited BED content
            provided as starting material.

        ::

            +--------------+-----------+-----------+--------------------+------------+
            | Chromosome   | Start     | End       | MeanCoverage       | Name       |
            | (category)   | (int32)   | (int32)   | (float64)          | (object)   |
            |--------------+-----------+-----------+--------------------+------------|
            | 1            | 155179779 | 155195266 | 274.39181248789305 | MUC1       |
            | 4            | 3072436   | 3079444   | 1060.255565068493  | HTT        |
            | 6            | 170554805 | 170563019 | 926.0217920623326  | SCA17      |
            | 9            | 27572666  | 27576436  | 521.1355437665783  | C9orf72    |
            | ...          | ...       | ...       | ...                | ...        |
            | X            | 147910539 | 147931116 | 597.4970112261262  | FMR1       |
            | 22           | 45791084  | 45797924  | 1638.6897660818713 | SCA10      |
            | 19           | 13205097  | 13210441  | 810.1785179640718  | SCA6       |
            | 14           | 92067404  | 92073500  | 1261.893372703412  | SCA3       |
            +--------------+-----------+-----------+--------------------+------------+

        >>> import megrim
        from megrim.environment import get_packaged_file_path

        bam = get_packaged_file_path()

        """
        logging.info("loading on target coordinates")
        on_target_universe = self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_target_ranges())
        on_target_universe.Name = self.bed.get_target_ranges().Name
        return on_target_universe

    @functools.lru_cache()
    def get_on_target_annotation(self):
        """
        Get pyranges format on-target coordinates augmented with bam info.

        This method returns annotation for the on_target
        regions of the genome. see
        :meth:`~megrim.target_enrichment.TargetEnrichment.get_on_target`
        The core coordinates are augmented with additional annotation scraped
        from the provided BAM file.

        Returns
        -------
        augmented_annotation: ``pyranges``
            A rich pyranges object augmented with additional data - see
            :meth:`~megrim.reference_genome.augment_annotation`

        """
        logging.info("loading on target annotation")
        augmented_annotation = augment_annotation(
            self.bam, self.get_on_target())
        return augmented_annotation

    @functools.lru_cache()
    def get_target_proximal(self):
        """
        Get pyranges format target-proximal coordinates and coverage.

        This method returns the coordinates for the target-proximal regions
        of the genome. Target proximity is defined at class initialisation
        with the ``target_proximity`` parameter. The targets are as
        defined with the provided bed file and as reported by
        :meth:`~megrim.target_enrichment.TargetEnrichment.get_on_target`

        Returns
        -------
        target_proximal_universe: ``pyranges``
            This is a pyranges refection of the regions of the genome that
            are immediately proximal to on-target regions.


        """
        logging.info("loading target proximal coordinates")
        return self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_target_proximal_ranges())

    @functools.lru_cache()
    def get_target_proximal_annotation(self):
        """
        Get contexually annotated info relating to target-proximal regions.

        See :meth:`~megrim.reference_genome.augment_annotation` for
        description of the method for extracting target proximal coordinates
        and basic metadata. This method takes these results and returns
        deeper context derived from a quick parse of the accompanying BAM
        file

        Returns
        -------
        augmented_annotation: ``pyranges``
            A rich pyranges object augmented with additional data - see
            :meth:`~megrim.reference_genome.augment_annotation`

        """
        logging.info("loading target proximal annotation")
        return augment_annotation(
            self.bam, self.get_target_proximal())

    @functools.lru_cache()
    def get_off_target(self, tile_size=None):
        """
        Get pyranges format off_target coordinates and coverage.

        This method returns the coordinates for the off-target but read-
        dense regions of the genome. This method builds a coverage-map of the
        genome and filters out the on-target and target proximal regions -
        the remaining untargetted genome is then assessed for mean coverage;
        the regions that correspond to >= mean.untargetted coverage *
        background_threshold are defined as being off target.

        Parameters
        ----------
        tile_size: ``int``
            This parameter defines the window size of the tiles to be
            placed across the genome. The smaller the window the better the
            resolution of off-target regions, but the greater the cost in
            compute time.

        Returns
        -------
        off_target_universe: ``pyranges``
            This is a pyranges refection of the regions of the genome that
            are defined as being off-target.
        """
        if tile_size is None:
            tile_size = self.tile_size
        logging.info(
            "loading off target coordinates (tile_size={})".format(tile_size))
        filtered_coverage = self.ref.get_tiled_mean_coverage(
            self.bam, tile_size=tile_size).subtract(
                self.get_on_target()).subtract(self.get_target_proximal())
        return self.ref.get_tiled_mean_coverage(
            self.bam, 
            ranges=filtered_coverage[
                filtered_coverage.MeanCoverage >= \
                    self.background_threshold * self.get_background_coverage()].slack(10).merge().slack(-10))
    
    @functools.lru_cache()
    def get_background_coverage(self, tile_size=None):
        """
        

        Parameters
        ----------
        tile_size : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        mean_coverage : TYPE
            DESCRIPTION.

        """
        if tile_size is None:
            tile_size = self.tile_size
        logging.info(
            "calculating background coverage (tile_size={})".format(tile_size))
        filtered_coverage = self.ref.get_tiled_mean_coverage(
            self.bam, tile_size=tile_size).subtract(
                self.get_on_target()).subtract(self.get_target_proximal())
        # this is an approximation since not all windows will be the same size
        mean_coverage = filtered_coverage.MeanCoverage.mean()
        return mean_coverage
            
    @functools.lru_cache()    
    def get_off_target_annotation(self, tile_size=None):    
        """
        

        Parameters
        ----------
        tile_size : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        logging.info("loading off target annotation")
        return augment_annotation(
            self.bam, self.get_off_target(tile_size=tile_size))
        
    @functools.lru_cache()
    def get_background(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        logging.info("loading background coordinates")
        return self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_untargeted_ranges().subtract(self.get_off_target()))
    
    @functools.lru_cache()
    def get_background_annotation(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        logging.info("loading background annotation")
        return augment_annotation(
            self.bam, self.get_background())
        
    def get_fine_coverage_aggregate(self, universe, target_proximity=None):
        """
        

        Parameters
        ----------
        universe : TYPE
            DESCRIPTION.
        target_proximity : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        aggregated_cov : TYPE
            DESCRIPTION.

        """
        if target_proximity is None:
            target_proximity = self.target_proximity
        aggregated_cov = self.ref.deep_dive(
            self.bam, universe, target_proximity=target_proximity)
        return aggregated_cov
    
    def get_mapped_bases(self): 
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.get_off_target_annotation().df.bases_start.sum() + \
            self.get_background_annotation().df.bases_start.sum() + \
            self.get_target_proximal_annotation().df.bases_start.sum() + \
            self.get_on_target_annotation().df.bases_start.sum()
            
    def get_on_target_perc(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.get_on_target_annotation().df.rstart.sum() / \
            (self.get_on_target_annotation().df.rstart.sum() + \
             self.get_off_target_annotation().df.rstart.sum() + \
                 self.get_background_annotation().df.rstart.sum() + \
                     self.get_target_proximal_annotation().df.rstart.sum()) * 100
                
    def get_depletion_factor(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        # depletion_label = "FAILED"
        return weighted_percentile(self.get_on_target_annotation()) / \
            weighted_percentile(self.get_background_annotation())
            
    def get_mean_target_coverage(self):
        """
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.get_on_target_annotation().df.MeanCoverage.repeat(
            self.get_on_target_annotation().End - \
                self.get_on_target_annotation().Start).mean()
            
    def get_cas9_exec_infographic(self, **kwargs):
        """
        

        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        (plot_width, plot_dpi) = self.handle_kwargs(
            ["plot_width", "plot_dpi"], **kwargs)
                
        throughput = InfographicNode(legend="Throughput",
                                        value="{:.2f} Gb".format(self.get_mapped_bases()/1e9),
                                        graphic='calculator')
        
        on_target_node = InfographicNode(legend="Reads on Target",
                                         value="{:.2f}%".format(self.get_on_target_perc()),
                                         graphic='cut')
        
        coverage_node = InfographicNode(legend="Mean Target Coverage",
                                     value="{:.2f}X".format(self.get_mean_target_coverage()),
                                     graphic='map-marked')
        
        depletion_node = InfographicNode(legend="Non-target Depletion",
                                     value="{:.1f}X".format(self.get_depletion_factor()),
                                     graphic='fill-drip')
        
        infographic_data = [throughput, on_target_node, coverage_node, depletion_node]
        ip = InfographicPlot(infographic_data, rows=1, columns=4)
        ip.plot_infographic(plot_width, plot_dpi)

        
    def get_mapping_stats(self, annotated_ranges, scale="Megabases"):
        """
        

        Parameters
        ----------
        annotated_ranges : TYPE
            DESCRIPTION.
        scale : TYPE, optional
            DESCRIPTION. The default is "Megabases".

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        scaleVal = 1
        if scale == "Gigabases":
            scaleVal = 1e9
        elif scale == "Megabases":
            scaleVal = 1e6
        elif scale == "Kilobases":
            scaleVal = 1e3
        
        total_reads = annotated_ranges.df.rstart.sum()
        fastq_nts = annotated_ranges.df.bases_start.sum() / scaleVal
        mapped_clipped = annotated_ranges.df.cigar_m.sum() / scaleVal
        genome_size = self.ref.get_genome_size()
        mapped_space = np.sum(annotated_ranges.End - annotated_ranges.Start)
        ref_space = mapped_space / genome_size * 100
        mean_cov = np.sum((annotated_ranges.End - annotated_ranges.Start) * \
            annotated_ranges.MeanCoverage) / \
            np.sum(annotated_ranges.End - annotated_ranges.Start)
        # return a pd.Series of observations
        return pd.Series({"Total reads": "{:.0f}".format(total_reads),
                          scale: "{:.2f}".format(fastq_nts),
                          "Mapped {}".format(scale): "{:.2f}".format(mapped_clipped),
                          "Genome fraction (%)": "{:.3f}".format(ref_space),
                          "Mean coverage (X)": "{:.2f}".format(mean_cov)})
    
    def get_mapping_summary_stats(self):
        """
        

        Returns
        -------
        df : TYPE
            DESCRIPTION.

        """
        on_t = self.get_mapping_stats(self.get_on_target_annotation())
        off_t = self.get_mapping_stats(self.get_off_target_annotation())
        t_p = self.get_mapping_stats(self.get_target_proximal_annotation())
        bg = self.get_mapping_stats(self.get_background_annotation())
        
        df = pd.concat([on_t, t_p, off_t, bg], axis=1, keys=["on_target", "target_proximal", "off_target", "background"])
        df.iloc[0] = df.iloc[0].astype(int)
        return df
            
    
    def get_target_performance(self, target, scale="Megabases"):
        """
        

        Parameters
        ----------
        target : TYPE
            DESCRIPTION.
        scale : TYPE, optional
            DESCRIPTION. The default is "Megabases".

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        scaleVal = 1
        if scale == "Gigabases":
            scaleVal = 1e9
        elif scale == "Megabases":
            scaleVal = 1e6
        elif scale == "Kilobases":
            scaleVal = 1e3
        df = self.get_on_target_annotation().df
        df = df.loc[df.Name==target,]
        
        return target, \
            "{:,}".format(int(int(df.End) - int(df.Start))), \
            "{:.2f}".format(float(df.MeanCoverage)), \
            "{:.0f}".format(float(df.rstart)), \
            "{:.2f}".format(float(df.bases_start) / scaleVal), \
            "{:.1f}".format(float(df.start_read_len)), \
            "{:.2f}".format(float(df.readq)), \
            "{:.2f}".format(float(df.mapq)), \
            "{:.1f}".format(float(df.strand_p / (df.strand_n + df.strand_p) * 100))    
        
    def get_target_performance_summary(self, scale="Megabases"):
        """
        

        Parameters
        ----------
        scale : TYPE, optional
            DESCRIPTION. The default is "Megabases".

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        targets = self.get_on_target().Name
        return targets.apply(self.get_target_performance).apply(
            pd.Series, 
            index=["Target", "Target size", "Mean coverage", "Read count", 
                   scale, "MeanReadLen", "MeanReadQ", "MeanMapQ", "(+)Strand"])
        
    def get_target_plot(self, target, **kwargs):
        """
        

        Parameters
        ----------
        target : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        targets = self.get_on_target().df
        targets = targets.loc[targets.Name == target,]
        
        aggregated_cov = self.ref.stranded_dive(
            bam=self.bam,
            arange=self.get_on_target(),
            target=target)
        #print(aggregated_cov)
        plot = figure(
            title="Plot showing depth of coverage across ({}) target region".format(targets.Name.tolist()[0]), 
            x_axis_label="Position on Chr({})".format(str(targets.Chromosome.tolist()[0])),
            y_axis_label='Depth of coverage (X)', 
            background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)
        
        plot.line(aggregated_cov.Start, aggregated_cov.MeanCoverage, 
                  line_width=2, line_color='#1F78B4', 
                  legend_label='Depth of Coverage')
        
        start_line = Span(location=targets.Start.tolist()[0], dimension='height', line_color='red', line_width=2, line_alpha=0.5)
        end_line = Span(location=targets.End.tolist()[0], dimension='height', line_color='red', line_width=2, line_alpha=0.5)
        bg_line = Span(location=(self.get_background_coverage()*self.background_threshold), dimension='width', line_color='orange', line_width=2, line_alpha=0.7)
        
        plot.renderers.extend([start_line, end_line, bg_line])
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)
        
    def get_stranded_plot(self, target, **kwargs):

        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        targets = self.get_on_target().df
        targets = targets.loc[targets.Name == target,]
        
        aggregated_cov = self.ref.stranded_dive(
            bam=self.bam,
            arange=self.get_on_target(),
            target=target).df
        aggregated_cov['baseline'] = 0
        aggregated_cov['topline'] = aggregated_cov['+'] + aggregated_cov['-']
        #print(aggregated_cov)
        plot = figure(
            title="Plot showing depth of coverage across ({}) target region".format(targets.Name.tolist()[0]), 
            x_axis_label="Position on Chr({})".format(str(targets.Chromosome.tolist()[0])),
            y_axis_label='Depth of coverage (X)', 
            background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)
        
        plot.varea_stack(stackers=["-","+"], x='Start', 
                         color=['#1F78B4', '#A6CEE3'], 
                         legend_label=["Forward strand", "Reverse strand"], 
                         source=aggregated_cov)
        
        start_line = Span(location=targets.Start.tolist()[0], dimension='height', line_color='red', line_width=2, line_alpha=0.5)
        end_line = Span(location=targets.End.tolist()[0], dimension='height', line_color='red', line_width=2, line_alpha=0.5)
        bg_line = Span(location=(self.get_background_coverage()*self.background_threshold), dimension='width', line_color='orange', line_width=2, line_alpha=0.5)
        
        plot.renderers.extend([start_line, end_line, bg_line])
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)
    
    def get_ideogram(self, **kwargs):
        """
        

        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        
        logging.info("preparing ideogram")
        
        chromosomes = self.ref.get_chromosomes()
        lengths = self.ref.get_chromosome_lengths(chromosomes)
        
        ideo_data = pd.DataFrame({"Chromosome": chromosomes,
                                  "Size": lengths}).reset_index(drop=True)
        # print(ideo_data)
        
        plot = figure(title='Ideogram showing location of off-target sequences', 
                      x_axis_label='Chromosome position (nt)',
                      y_axis_label='Chromosome', background_fill_color="lightgrey",
                      plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        plot.hbar(y=ideo_data.index.tolist(), right=ideo_data.Size, left=0, height=0.7, fill_color="#A6CEE3", line_color="#1F78B4")
        
        # and overlay some coordinate data
        off_t = self.get_off_target().df
        def getit(i):
            return ideo_data.loc[ideo_data["Chromosome"]==i,].index.tolist()[0]
        
        off_t['row_id'] = off_t['Chromosome'].apply(getit)
        
        plot.hbar(y=off_t['row_id'], left=off_t['Start'], right=off_t['End'], height=0.7, fill_color="red", line_color="red")
        
        plot.yaxis.ticker = ideo_data.index.tolist()
        tick_dict = {}
        for i in ideo_data.index.tolist():
            tick_dict[i]=ideo_data.iloc[i].Chromosome
        plot.yaxis.major_label_overrides = tick_dict
        #return ideo_data
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)
        # return ideo_data
    
    
    def get_off_target_stats(self):
        """
        

        Returns
        -------
        df : TYPE
            DESCRIPTION.

        """
        data = self.get_off_target_annotation().df
    
        df = pd.concat(
            [data["Chromosome"],
             data["Start"].apply("{:,}".format),
             data["End"].apply("{:,}".format),
             (data["End"] - data["Start"]).apply("{:,}".format),
             data["MeanCoverage"].apply("{:.2f}".format),
             (data.strand_n + data.strand_p),
             data["mean_read_len"].apply("{:.2f}".format),
             (data.strand_p / (data.strand_n + data.strand_p) * 100).apply("{:.2f}".format),
             data["read0"].apply("{:.2f}".format),
             data["map0"].apply("{:.2f}".format)
             ], axis=1, 
            keys=["Chromosome", "Start", "End", "Width", "MeanCoverage", 
                  "MappedReads", "MeanReadLength", "%FWD", "MeanReadQ", 
                  "MeanMapQ"])
 
        return df
      
    
    