#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:30:22 2020.

@author: srudd
"""

from megrim.environment import Flounder
from megrim.reference_genome import ReferenceGenome
from megrim.genome_geometry import BamHandler
from bokeh.plotting import figure
import numpy as np
import pandas as pd
import pyranges as pr
from bokeh.models import NumeralTickFormatter, Label


class VirusGenome(Flounder):
    """
    Class for methods relating to plotting depth of coverage for small genomes.

    This is a collection of methods to accessorise a workflow looking at
    viral genomes. The code in this workflow does not really need to stay
    here and should probably be merged into the ReferenceGenome class.
    """

    # todo: please merge this class into the ReferenceGenome class.

    def __init__(self, ref, bam, bam_b=None, fasta=None):
        Flounder.__init__(self)
        self.ref = ReferenceGenome(ref)
        self.bam = BamHandler(bam)
        self.bam_b = None
        if bam_b is not None:
            self.bam_b = BamHandler(bam_b)
        self.fasta = fasta

    def get_coverage(self, plot_bam_b=False, tile_size=10):
        """
        Get the mean depth-of-coverage across the specified genome.

        The
        :meth:`~megrim.reference_genome.ReferenceGenome.get_tiled_mean_coverage`
        framework is used in combination with the depth-of-coverage
        information in the BAM file to look for depth of coverage across the
        genome using the specified tile size. This method simply returns the
        basic coverage information.

        Parameters
        ----------
        plot_bam_b: boolean
            This specifies whether a supplementary RG bam should be plotted instead.
            The default is False (standard bam).
        tile_size: int
            The size of the window to use when tiling the genome.
            The default is 10.

        Returns
        -------
        pandas.DataFrame
            A DataFrame describing the bins across the genome, their
            boundaries and mean depths-of-coverage.

        """
        if plot_bam_b and self.bam_b is not None:
            return self.ref.get_tiled_mean_coverage(
                self.bam_b, tile_size=tile_size)
        return self.ref.get_tiled_mean_coverage(
            self.bam, tile_size=tile_size)

    def plot_coverage(self, tile_size=10, id="undefined", **kwargs):
        """
        Prepare a coverage plot across a whole chromosome.

        This method is assuming that we are working with a single chromosome.
        The get_coverage method is used to prepare tiles across the genome
        and they are displayed in a plot. The plot dimensions and output are
        configurable through the **kwargs.

        Parameters
        ----------
        tile_size: int
            The size of the window to use when tiling the genome.
            The default is 10.
        id: String
            To be used in the presentation of the figure title
        **kwargs: **kwargs
            can provide a number of possible options to the Flounder class in
            the background that may be used to alter the plot dimensions,
            bokeh tools and plot rendering options.

        Returns
        -------
        bokeh image plot
            The plot returned to the calling method - please check
            :meth:`~megrim.environment.Flounder` for more information on
            how this can be configured.

        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)
        coverage = self.get_coverage(tile_size=tile_size)
        coverage_b = None
        if self.bam_b is not None:
            coverage_b = self.get_coverage(plot_bam_b=True, tile_size=tile_size)

        plot = figure(
            title=f"Plot showing depth of coverage across ({id}) genome",
            x_axis_label="Position on Chr({})".
            format(str(coverage.Chromosome.tolist()[0])),
            y_axis_label='Depth of coverage (X)',
            background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        if coverage_b is not None:
            plot.line(coverage.Start, coverage.MeanCoverage, line_width=2, line_color='#1F78B4',
                      legend_label='Depth-of-Coverage (RG1)')
            plot.line(coverage_b.Start, coverage_b.MeanCoverage, line_width=2, line_color='#A6CEE3',
                      legend_label='Depth-of-Coverage (RG2)')
        else:
            plot.line(coverage.Start, coverage.MeanCoverage, line_width=2, line_color='#1F78B4',
                      legend_label='Depth-of-Coverage')

        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        return self.handle_output(plot, plot_type)

    def plot_coverage_distribution(
            self, tile_size=10, bins=30, id="undefined", **kwargs):
        """
        Plot a histogram showing the distribution of depths-of-coverage.

        This plot breaks the genome into bins of tile_size in size. The
        depth-of-coverage across these bins is then binned to yield a
        histogram showing the relative frequency of different depths of
        coverage across the genome.

        Parameters
        ----------
        tile_size: int
            The size of the window to use when tiling the genome.
            The default is 10.
        bins: int
            The number of bins to break the distribution into.
            The default is 30.
        id: String
            To be used in the presentation of the figure title
        **kwargs: **kwargs
            can provide a number of possible options to the Flounder class in
            the background that may be used to alter the plot dimensions,
            bokeh tools and plot rendering options.

        Returns
        -------
        bokeh image plot
            The plot returned to the calling method - please check
            :meth:`~megrim.environment.Flounder` for more information on
            how this can be configured.
        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)

        coverage = self.get_coverage(tile_size=tile_size)
        if self.bam_b is not None:
            coverage = pd.merge(self.get_coverage(tile_size=tile_size).df,
                                self.get_coverage(tile_size=tile_size, plot_bam_b=True).df,
                                on=["Chromosome", "Start", "End"])
            coverage['MeanCoverage'] = coverage[["MeanCoverage_x", "MeanCoverage_y"]].sum(axis=1)

        mc = coverage.MeanCoverage
        deepest_bin = mc.max()+1
        boundaries = np.linspace(
            0, deepest_bin, num=bins, endpoint=True, retstep=False)
        assignments = np.digitize(mc, boundaries)

        xxx = pd.DataFrame({"coverage": mc,
                            "assignment": assignments,
                            "tally": tile_size}).groupby(
                                ["assignment"]).agg(
                                    {"assignment": ["first"],
                                     "tally": [np.sum]})
        xxx.columns = xxx.columns.droplevel()
        yyy = pd.DataFrame({"start": boundaries[:-1],
                            "end": boundaries[1:]},
                           index=np.arange(1, len(boundaries)))
        coverage_dist = yyy.merge(
            xxx, how="outer", left_index=True, right_index=True)
        coverage_dist.columns = ["start", "end", "batch", "count"]
        coverage_dist["batch"] = coverage_dist.index
        coverage_dist = coverage_dist.fillna(0)
        coverage_dist["colour"] = "#1F78B4"

        # print(coverage_dist)
        # print(coverage_dist['count'].sum())

        p = figure(title="Histogram showing distribution of coverage",
                   background_fill_color="lightgrey", plot_width=plot_width,
                   plot_height=plot_height, tools=plot_tools)
        p.quad(
            source=coverage_dist, top="count", bottom=0, left='start',
            right='end', fill_color='colour', line_color="white", alpha=0.7)
        p.xaxis.axis_label = 'Depth-of-coverage (X-fold)'
        p.yaxis.axis_label = 'Bases of genome (n)'
        return self.handle_output(p, plot_type)

    def plot_coverage_fractions(self, max_depth=None, bins=100, tile_size=10,
                                boundaries_of_interest=[10, 100], **kwargs):
        """
        Prepare ROC-like plot showing distribution of coverage across genome.

        This method prepares a plot that aims to visualise the fraction of
        the genome that is reflected at different levels of coverage. The
        genome is split into windows of size tile_size and coverage averaged
        across these windows. The maximum depth is calculated and the
        coverages are binned into bins bins. The data is then prepared so
        that 100% of the genome is covered at 0 depth going to 0% of the
        genome at max_depth+1 depth.

        Parameters
        ----------
        max_depth: int
            The maximum depth to render within the plot.
            # TODO: this looks borked
        bins: int
            The number of bins to present in the plot. This is really an
            aesthetic thing only.
        tile_size: int
            The tile size to summarise depths at. The default is 10.
        boundaries_of_interest: list of ints
            This will direct the plotting of lines to show the positions where
            specific coverages are reached. The default is [10, 100].
        **kwargs: **kwargs
            can provide a number of possible options to the Flounder class in
            the background that may be used to alter the plot dimensions,
            bokeh tools and plot rendering options.

        Returns
        -------
        bokeh image plot

        """
        (plot_width, plot_height, plot_type, plot_tools) = self.handle_kwargs(
            ["plot_width", "plot_height", "plot_type", "plot_tools"], **kwargs)

        coverage = self.get_coverage(tile_size=tile_size)
        if self.bam_b is not None:
            coverage = pd.merge(self.get_coverage(tile_size=tile_size).df,
                                self.get_coverage(tile_size=tile_size, plot_bam_b=True).df,
                                on=["Chromosome", "Start", "End"])
            coverage['MeanCoverage'] = coverage[["MeanCoverage_x", "MeanCoverage_y"]].sum(axis=1)
            # and convert back to pyranges ...
            coverage = pr.PyRanges(coverage)

        if max_depth is None:
            max_depth = coverage.MeanCoverage.max() + 1
            print("max_depth set to {}".format(max_depth))
        boundaries = np.linspace(
            0, max_depth, num=bins, endpoint=True, retstep=False)
        assignments = np.digitize(coverage.MeanCoverage, boundaries)
        cov_data = pd.DataFrame({"assignment": assignments,
                                 "bases": (coverage.End - coverage.Start)})
        cov_data = cov_data.groupby("assignment").agg({"assignment": "first",
                                                       "bases": np.sum})
        # add the missing values
        cov_data = cov_data.reindex(
            pd.Index(
                pd.Series(boundaries).index)
            ).reset_index().drop(["assignment"], axis=1).fillna(0)
        # add the cumulative sum to the data
        cov_data['frac'] = cov_data.bases.sum() - np.cumsum(cov_data.bases)
        # prepare the cumsum as percentage
        cov_data['perc'] = cov_data.frac / cov_data.bases.sum() * 100

        plot = figure(
            title='Plot showing % of genome covered at different depths',
            x_axis_label='Depth of coverage (X)',
            y_axis_label='Percentage of Genome covered(%)',
            background_fill_color="lightgrey",
            plot_width=plot_width, plot_height=plot_height, tools=plot_tools)

        for b in boundaries_of_interest:
            # b is coverage - get corresponding perc
            bases = coverage[coverage.MeanCoverage >= b].lengths().sum()
            perc = bases / cov_data.bases.sum() * 100
            legend = "{}X".format(b)
            plot.line(
                [0, b, b], [perc, perc, 0], line_width=2, line_color='red')
            plot.add_layout(Label(x=b, y=perc, text=legend, text_color='red'))

        plot.step(boundaries, cov_data.perc, line_width=2, mode="before")
        return self.handle_output(plot, plot_type)
