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

        Returns:
        -------
        on_target_universe: ``pyranges``
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

        See
        :meth:`~megrim.target_enrichment.TargetEnrichent.get_target_proximal`
        for description of method for extracting target proximal coordinates
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
                filtered_coverage.MeanCoverage >=
                self.background_threshold *
                self.get_background_coverage()].slack(
                    10).merge().slack(-10))

    @functools.lru_cache()
    def get_background_coverage(self, tile_size=None):
        """
        Get mean background depth of coverage across the genome.

        Background refers to the part of the genome that is not on-target,
        target-proximal or off-target. This method tiles the reference
        genome and removes the on-target and target-proximal regions to
        identify the background coverage. The coverage across these remaining
        regions is averaged to return the background. This method can be
        tuned by tile_size that defines the size of the tiles to place
        across the genome; smaller tiles better resolution but computationally
        more demanding

        .. note::

            Background is calculated using a heuristic that assumes that the
            window size is equal in all cases ... this approximation may
            not be ideal ... please shout if this should be fixed!

        Parameters
        ----------
        tile_size: INT
            The size of the tiles to place across the genome.

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
        Get contexually annotated info relating to off_target regions.

        See :meth:`~megrim.target_enrichment.TargetEnrichent.get_off_target`
        for description of method for extracting off_target coordinates
        and basic metadata. This method takes these results and returns
        deeper context derived from a quick parse of the accompanying BAM
        file

        Parameters
        ----------
        tile_size: ``int``
            This parameter defines the window size of the tiles to be
            placed across the genome. The smaller the window the better the
            resolution of off-target regions, but the greater the cost in
            compute time.

        Returns
        -------
        augmented_annotation: ``pyranges``
            A rich pyranges object augmented with additional data - see
            :meth:`~megrim.reference_genome.augment_annotation`
        """
        logging.info("loading off target annotation")
        return augment_annotation(
            self.bam, self.get_off_target(tile_size=tile_size))

    @functools.lru_cache()
    def get_background(self):
        """
        Get pyranges format background coordinates and coverage.

        Background refers to the part of the genome that is not on-target,
        target-proximal or off-target. This method tiles the reference
        genome and removes the on-target and target-proximal regions to
        calculate background coverage. The off-target regions are identified
        and subtracted to yield the background.

        Returns
        -------
        background_universe: ``pyranges``
            This is a pyranges refection of the regions of the genome that
            are defined as being background.

        """
        logging.info("loading background coordinates")
        return self.ref.get_tiled_mean_coverage(
            self.bam, ranges=self.bed.get_untargeted_ranges().subtract(
                self.get_off_target()))

    @functools.lru_cache()
    def get_background_annotation(self):
        """
        Get contexually annotated info relating to background regions.

        See :meth:`~megrim.target_enrichment.TargetEnrichent.get_background`
        for description of method for extracting background coordinates
        and basic metadata. This method takes these results and returns
        deeper context derived from a quick parse of the accompanying BAM
        file

        Parameters
        ----------
        tile_size: ``int``
            This parameter defines the window size of the tiles to be
            placed across the genome. The smaller the window the better the
            resolution of off-target regions, but the greater the cost in
            compute time.

        Returns
        -------
        augmented_annotation: ``pyranges``
            A rich pyranges object augmented with additional data - see
            :meth:`~megrim.reference_genome.augment_annotation`
        """
        logging.info("loading background annotation")
        return augment_annotation(
            self.bam, self.get_background())

    def get_fine_coverage_aggregate(self, universe,
                                    window_size=10, target_proximity=None):
        """
        Prepare finer resolution coverage information for regions of interest.

        The methods such as
        :meth:`~megrim.target_enrichment.TargetEnrichent.get_on_target_annotation`
        report sequence characteristics for genomic ranges of interest;
        key values are summarised. For the preparation of graphs and other
        analyses a finer resolution view of a region of interest should be
        performed - this method takes existing ranges and re-windows with a
        smaller interval (5nt, 10nt ...) as defined by window_size.

        Parameters
        ----------
        universe: ``pyranges``
            The starting coordinates to further refine.
        window_size: ``int``
            The size of the window to apply to the provided ranges
        target_proximity : ``int``
            The regions immediately up- and down-stream of region of interest

        Returns
        -------
        aggregated_cov: ``pyranges``
            A finer resolution pyranges object tiled over the provided
            ranges
        """
        if target_proximity is None:
            target_proximity = self.target_proximity
        aggregated_cov = self.ref.deep_dive(
            self.bam, universe, target_proximity=target_proximity)
        return aggregated_cov

    def get_mapped_bases(self):
        """
        Return number of mapped based described within the contained ranges.

        This is an accessory method to report the number of sequence bases
        that are contained within the background, off-target, target-proximal
        and on-target regions of the genome.

        Returns
        -------
        mapped_bases: ``int``
            The sum of mapped bases across the universes considered.

        """
        return self.get_off_target_annotation().df.bases_start.sum() + \
            self.get_background_annotation().df.bases_start.sum() + \
            self.get_target_proximal_annotation().df.bases_start.sum() + \
            self.get_on_target_annotation().df.bases_start.sum()

    def get_on_target_perc(self):
        """
        Return the percentage of sequence reads that correspond to on-target.

        This accessory method returns the percentage of sequence reads that
        map to on-target ranges. An on-target read is here defined as a read
        that either starts or ends within the given target range.

        Returns
        -------
        on_target_perc: ``float``
            Percentage of on-target reads.
        """
        return self.get_on_target_annotation().df.rstart.sum() / \
            (self.get_on_target_annotation().df.rstart.sum() +
             self.get_off_target_annotation().df.rstart.sum() +
             self.get_background_annotation().df.rstart.sum() +
             self.get_target_proximal_annotation().df.rstart.sum()) * 100

    def get_depletion_factor(self):
        """
        Calculate the depletion factor following target-enrichment protocol.

        A target enrichment may in fact deplete the non-target regions.
        Mind-blown! This method here returns a simple estimate of overall
        of non-target depletion by returning the 50 percentile of
        on-target depth-of-coverage / 50 percentile of background coverage.
        This means that off-target regions are excluded from this calculation.

        Returns
        -------
        depletion_factor: float
            The depletion factor.

        """
        # depletion_label = "FAILED"
        return weighted_percentile(self.get_on_target_annotation()) / \
            weighted_percentile(self.get_background_annotation())

    def get_mean_target_coverage(self):
        """
        Return mean sequence coverage over target regions.

        This accessory method calculates the mean coverage across the
        on-target regions defined by the method
        :meth:`~megrim.target_enrichment.TargetEnrichment.get_on_target`

        Returns
        -------
        mean_coverage: ``int``
            The mean coverage across all provided target regions.

        """
        return self.get_on_target_annotation().df.MeanCoverage.repeat(
            self.get_on_target_annotation().End -
            self.get_on_target_annotation().Start).mean()

    def get_cas9_exec_infographic(self, **kwargs):
        """
        Prepare infographic-like plot summarising key enrichment metrics.

        This function returns an infographic-like plot that summaries the
        performance of the target enrichment analysis.

        Parameters
        ----------
        **kwargs: **kwargs
            can provide a number of possible options to the Flounder class in
            the background that may be used to alter the plot_scaling and
            dpi values.

        Returns
        -------
        None.
            A plot will be displayed on the current graphics device - this
            should probably be reviewed for greater consistency with the
            bokeh styled plots

        """
        (plot_width, plot_dpi) = self.handle_kwargs(
            ["plot_width", "plot_dpi"], **kwargs)

        throughput = InfographicNode(
            legend="Throughput",
            value="{:.2f} Gb".format(self.get_mapped_bases()/1e9),
            graphic='calculator')

        on_target_node = InfographicNode(
            legend="Reads on Target",
            value="{:.2f}%".format(self.get_on_target_perc()),
            graphic='cut')

        coverage_node = InfographicNode(
            legend="Mean Target Coverage",
            value="{:.2f}X".format(self.get_mean_target_coverage()),
            graphic='map-marked')

        depletion_node = InfographicNode(
            legend="Non-target Depletion",
            value="{:.1f}X".format(self.get_depletion_factor()),
            graphic='fill-drip')

        infographic_data = [
            throughput, on_target_node, coverage_node, depletion_node]
        ip = InfographicPlot(infographic_data, rows=1, columns=4)
        ip.plot_infographic(plot_width, plot_dpi)

    def get_mapping_stats(self, annotated_ranges, scale="Megabases"):
        """
        Prepare basic mapping statistics from a provided pyranges.

        This accessory method prepares basic mapping statistics for the
        provided pyranges object - this can be used to summarise mapping
        for given regions of interest.

        Parameters
        ----------
        annotated_ranges: ``pyranges``
            This is the range of interest - could be e.g. on_target from
            :meth:`~megrim.target_enrichment.TargetEnrichent.get_on_target`

        scale: ``str``
            The default is "Megabases". Also allowed are ["Gigabases",
            "Kilobases"]

        Returns
        -------
        mapping_stats: ``pandas.Series``
            try me - this is used for functions such as
            :meth:`~megrim.target_enrichment.TargetEnrichent.get_mapping_summary_stats`

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
        mean_cov = np.sum(((
            annotated_ranges.End - annotated_ranges.Start) *
            annotated_ranges.MeanCoverage) /
            np.sum(annotated_ranges.End - annotated_ranges.Start))
        # return a pd.Series of observations
        return pd.Series({
            "Total reads": "{:.0f}".format(total_reads),
            scale: "{:.2f}".format(fastq_nts),
            "Mapped {}".format(scale): "{:.2f}".format(mapped_clipped),
            "Genome fraction (%)": "{:.3f}".format(ref_space),
            "Mean coverage (X)": "{:.2f}".format(mean_cov)})

    def get_mapping_summary_stats(self):
        """
        Prepare mapping statistics across the whole dataset.

        This function prepares mapping statistics for display in e.g.
        Jupyter notebooks that summarise the mapping characteristics for the
        canonical target-universes considered in an enrichment study. The
         :meth:`~megrim.target_enrichment.TargetEnrichent.get_mapping_stats`
         function is called on the canonical levels.

        Returns
        -------
        mapping_summary_stats: ``Pandas.DataFrame``
            A DataFrame object summarising the key mapping observations.

        """
        on_t = self.get_mapping_stats(self.get_on_target_annotation())
        off_t = self.get_mapping_stats(self.get_off_target_annotation())
        t_p = self.get_mapping_stats(self.get_target_proximal_annotation())
        bg = self.get_mapping_stats(self.get_background_annotation())

        df = pd.concat(
            [on_t, t_p, off_t, bg],
            axis=1,
            keys=["on_target", "target_proximal", "off_target", "background"])
        df.iloc[0] = df.iloc[0].astype(int)
        return df

    def get_target_performance(self, target, scale="Megabases"):
        """
        Prepare mapping characteristics for a target region of interest.

        This function summarises the mapping performance for a single
        range corresponding to a named target of interest. Summary statistics
        such as the number of reads, number of bases and strandedness for the
        mapped reads is reported along with mean read quality and mapping
        quality.

        Parameters
        ----------
        target: ``Str``
            The target of interest to report.
        scale: ``Str``
            A variable to define how the bases should be summarised.
            The default is "Megabases". Also allowed are ["Gigabases",
            "Kilobases"]

        Returns
        -------
        tuple
            An unnamed tuple containing the elements of interest.

        """
        scaleVal = 1
        if scale == "Gigabases":
            scaleVal = 1e9
        elif scale == "Megabases":
            scaleVal = 1e6
        elif scale == "Kilobases":
            scaleVal = 1e3
        df = self.get_on_target_annotation().df
        df = df.loc[df.Name == target, ]

        return target, \
            "{:,}".format(int(int(df.End) - int(df.Start))), \
            "{:.2f}".format(float(df.MeanCoverage)), \
            "{:.0f}".format(float(df.rstart)), \
            "{:.2f}".format(float(df.bases_start) / scaleVal), \
            "{:.1f}".format(float(df.start_read_len)), \
            "{:.2f}".format(float(df.readq)), \
            "{:.2f}".format(float(df.mapq)), \
            "{:.1f}".format(float(
                df.strand_p / (df.strand_n + df.strand_p) * 100))

    def get_target_performance_summary(self, scale="Megabases"):
        """
        Prepare mapping characteristics for all target regions of interest.

        This accessory method collates information from the
        :meth:`~megrim.target_enrichment.TargetEnrichent.get_target_performance`
        function into a Pandas.DataFrame for display in downstream
        applications.

        Parameters
        ----------
        scale: ``Str``
            A variable to define how the bases should be summarised.
            The default is "Megabases". Also allowed are ["Gigabases",
            "Kilobases"]

        Returns
        -------
        mapping_summary_stats: ``Pandas.DataFrame``
            A DataFrame object summarising the key mapping observations.

        """
        targets = self.get_on_target().Name
        return targets.apply(self.get_target_performance).apply(
            pd.Series,
            index=["Target", "Target size", "Mean coverage", "Read count",
                   scale, "MeanReadLen", "MeanReadQ", "MeanMapQ", "(+)Strand"])

    def get_target_plot(self, target, **kwargs):
        """
        Prepare depth-of-coverage plot for a target region of interest.

        This function prepares a line plot that summarises the depth of
        coverage over a named target region of interest. The plot includes
        coverage over the target region expanded to include the up- and down-
        stream regions as enumerated with the target_proximity variable.

        Parameters
        ----------
        target: ``Str``
            This is the name of the target region and must correspond to a
            named object as seen within the results from
             :meth:`~megrim.target_enrichment.TargetEnrichment.get_on_target`.
        **kwargs: **kwargs
            can provide a number of possible options to the Flounder class in
            the background that may be used to alter the plot dimensions,
            bokeh tools and plot rendering options.

        Returns
        -------
        plot
            A bokeh derived plot that may be in a variery of different forms.

        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)

        targets = self.get_on_target().df
        targets = targets.loc[targets.Name == target, ]

        aggregated_cov = self.ref.stranded_dive(
            bam=self.bam,
            arange=self.get_on_target(),
            target=target)
        # print(aggregated_cov)
        plot = figure(
            title="Plot showing depth of coverage across ({}) target region".
            format(targets.Name.tolist()[0]),
            x_axis_label="Position on Chr({})".
            format(str(targets.Chromosome.tolist()[0])),
            y_axis_label='Depth of coverage (X)',
            background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        plot.line(aggregated_cov.Start, aggregated_cov.MeanCoverage,
                  line_width=2, line_color='#1F78B4',
                  legend_label='Depth of Coverage')

        start_line = Span(
            location=targets.Start.tolist()[0],
            dimension='height', line_color='red', line_width=2,
            line_alpha=0.5)
        end_line = Span(
            location=targets.End.tolist()[0], dimension='height',
            line_color='red', line_width=2, line_alpha=0.5)
        bg_line = Span(
            location=(
                self.get_background_coverage()*self.background_threshold),
            dimension='width', line_color='orange', line_width=2,
            line_alpha=0.7)

        plot.renderers.extend([start_line, end_line, bg_line])
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)

    def get_stranded_plot(self, target, **kwargs):
        """
        Prepare stranded depth-of-coverage plot for target of interest.

        This function prepares a stacked varea plot that summarises the
        depth of coverage for reads that map to the forward and reverse
        strands of the reference sequence for a named target region of
        interest. The plotred range includes coverage over the target region
        expanded to include the up- and down-stream regions as enumerated
        with the target_proximity variable.

        Parameters
        ----------
        target: ``Str``
            This is the name of the target region and must correspond to a
            named object as seen within the results from
             :meth:`~megrim.target_enrichment.TargetEnrichment.get_on_target`.
        **kwargs: **kwargs
            can provide a number of possible options to the Flounder class in
            the background that may be used to alter the plot dimensions,
            bokeh tools and plot rendering options.

        Returns
        -------
        plot
            A bokeh derived plot that may be in a variery of different forms.

        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)

        targets = self.get_on_target().df
        targets = targets.loc[targets.Name == target, ]

        aggregated_cov = self.ref.stranded_dive(
            bam=self.bam,
            arange=self.get_on_target(),
            target=target).df
        aggregated_cov['baseline'] = 0
        aggregated_cov['topline'] = aggregated_cov['+'] + aggregated_cov['-']
        # print(aggregated_cov)
        plot = figure(
            title="Plot showing depth of coverage across ({}) target region".
            format(targets.Name.tolist()[0]),
            x_axis_label="Position on Chr({})".
            format(str(targets.Chromosome.tolist()[0])),
            y_axis_label='Depth of coverage (X)',
            background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        plot.varea_stack(stackers=["-", "+"], x='Start',
                         color=['#1F78B4', '#A6CEE3'],
                         legend_label=["Forward strand", "Reverse strand"],
                         source=aggregated_cov)

        start_line = Span(
            location=targets.Start.tolist()[0], dimension='height',
            line_color='red', line_width=2, line_alpha=0.5)
        end_line = Span(
            location=targets.End.tolist()[0], dimension='height',
            line_color='red', line_width=2, line_alpha=0.5)
        bg_line = Span(
            location=(
                self.get_background_coverage()*self.background_threshold),
            dimension='width', line_color='orange', line_width=2,
            line_alpha=0.5)

        plot.renderers.extend([start_line, end_line, bg_line])
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)

    def get_ideogram(self, **kwargs):
        """
        Prepare an ideogram showing the genomic location of off-target regions.

        An ideogram summarises chromosomes ranked by name, sized by
        chromosomal length and overlaid with spatial information that in this
        case describes the genomic locations for regions of the genome that
        contain off-target peaks of sequence mapping.

        Parameters
        ----------
        **kwargs: **kwargs
            can provide a number of possible options to the Flounder class in
            the background that may be used to alter the plot dimensions,
            bokeh tools and plot rendering options.

        Returns
        -------
        plot
            A bokeh derived plot that may be in a variery of different forms.

        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)

        logging.info("preparing ideogram")

        chromosomes = self.ref.get_chromosomes()
        lengths = self.ref.get_chromosome_lengths(chromosomes)

        ideo_data = pd.DataFrame({"Chromosome": chromosomes,
                                  "Size": lengths}).reset_index(drop=True)
        # print(ideo_data)

        plot = figure(
            title='Ideogram showing location of off-target sequences',
            x_axis_label='Chromosome position (nt)',
            y_axis_label='Chromosome', background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        plot.hbar(
            y=ideo_data.index.tolist(), right=ideo_data.Size,
            left=0, height=0.7, fill_color="#A6CEE3", line_color="#1F78B4")

        # and overlay some coordinate data
        off_t = self.get_off_target().df

        def getit(i):
            return ideo_data.loc[
                ideo_data["Chromosome"] == i, ].index.tolist()[0]

        off_t['row_id'] = off_t['Chromosome'].apply(getit)

        plot.hbar(
            y=off_t['row_id'], left=off_t['Start'], right=off_t['End'],
            height=0.7, fill_color="red", line_color="red")

        plot.yaxis.ticker = ideo_data.index.tolist()
        tick_dict = {}

        for i in ideo_data.index.tolist():
            tick_dict[i] = ideo_data.iloc[i].Chromosome

        plot.yaxis.major_label_overrides = tick_dict
        # return ideo_data
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)
        # return ideo_data

    def get_off_target_stats(self):
        """
        Prepare tabular data summarising off-target regions of the genome.

        This method summarises some of the mapping characteristics for the
        off-target regions of the genome. As with
        :meth:`~megrim.target_enrichment.TargetEnrichent.get_target_performance_summary`
        the reported data includes information on mapping quality,
        read quality and strandedness of mapping.

        Returns
        -------
        df: ``Pandas.DataFrame``
            The results in tabular format for downstream export to workbooks
            and spreadsheet files.

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
             (data.strand_p / (
                 data.strand_n + data.strand_p) * 100).apply("{:.2f}".format),
             data["read0"].apply("{:.2f}".format),
             data["map0"].apply("{:.2f}".format)
             ], axis=1,
            keys=["Chromosome", "Start", "End", "Width", "MeanCoverage",
                  "MappedReads", "MeanReadLength", "%FWD", "MeanReadQ",
                  "MeanMapQ"])
        return df
