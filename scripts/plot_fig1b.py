import argparse
import copy
import errno
import io
import itertools
import logging
import os
import pickle
import re
import string
import sys
import warnings
from bdb import set_trace
from collections import OrderedDict


import matplotlib
import matplotlib.font_manager
import matplotlib.lines as lines
import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import sklearn
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LinearLocator
from PIL import Image
from sklearn.preprocessing import LabelEncoder

import sigProfilerPlotting as spplt

matplotlib.use("Agg")


# import sigProfilerPlotting as sigPlt
def rev_comp(mut: str):
    translator = str.maketrans("ACGT", "TGCA")
    mut = mut[-1] + mut[1:-1] + mut[0]
    new_mut = mut.translate(translator)
    return new_mut

from sigProfilerPlotting import (
    load_custom_fonts, process_input, reindex_sbs96, 
    SPP_TEMPLATES, getylabels, output_results,
)


def make_pickle_file_SBS96(context="SBS96", return_plot_template=False, volume=None, revcc=False):
    if revcc:
        print('revcc')
    if volume is None:
        volume = SPP_TEMPLATES

    path = os.path.join(volume, context + ".pkl")
    print(path)

    # if the pickle file already exists, return the template
    # if os.path.exists(path):
    #     return pickle.load(open(path, "rb"))

    # check if the template directory exists, create if not
    if not os.path.exists(volume):
        os.mkdir(volume)

    if context == "SBS96":
        plot_custom_text = False
        sig_probs = False
        pcawg = False

        # total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
        # , extent=[-5, 80, -5, 30])
        plt.rcParams["axes.linewidth"] = 2
        plot1 = plt.figure(figsize=(43.93, 9.92))
        plt.rc("axes", edgecolor="lightgray")
        panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
        seq96 = [
            "A[C>A]A",
            "A[C>A]C",
            "A[C>A]G",
            "A[C>A]T",
            "C[C>A]A",
            "C[C>A]C",
            "C[C>A]G",
            "C[C>A]T",
            "G[C>A]A",
            "G[C>A]C",
            "G[C>A]G",
            "G[C>A]T",
            "T[C>A]A",
            "T[C>A]C",
            "T[C>A]G",
            "T[C>A]T",
            "A[C>G]A",
            "A[C>G]C",
            "A[C>G]G",
            "A[C>G]T",
            "C[C>G]A",
            "C[C>G]C",
            "C[C>G]G",
            "C[C>G]T",
            "G[C>G]A",
            "G[C>G]C",
            "G[C>G]G",
            "G[C>G]T",
            "T[C>G]A",
            "T[C>G]C",
            "T[C>G]G",
            "T[C>G]T",
            "A[C>T]A",
            "A[C>T]C",
            "A[C>T]G",
            "A[C>T]T",
            "C[C>T]A",
            "C[C>T]C",
            "C[C>T]G",
            "C[C>T]T",
            "G[C>T]A",
            "G[C>T]C",
            "G[C>T]G",
            "G[C>T]T",
            "T[C>T]A",
            "T[C>T]C",
            "T[C>T]G",
            "T[C>T]T",
            "A[T>A]A",
            "A[T>A]C",
            "A[T>A]G",
            "A[T>A]T",
            "C[T>A]A",
            "C[T>A]C",
            "C[T>A]G",
            "C[T>A]T",
            "G[T>A]A",
            "G[T>A]C",
            "G[T>A]G",
            "G[T>A]T",
            "T[T>A]A",
            "T[T>A]C",
            "T[T>A]G",
            "T[T>A]T",
            "A[T>C]A",
            "A[T>C]C",
            "A[T>C]G",
            "A[T>C]T",
            "C[T>C]A",
            "C[T>C]C",
            "C[T>C]G",
            "C[T>C]T",
            "G[T>C]A",
            "G[T>C]C",
            "G[T>C]G",
            "G[T>C]T",
            "T[T>C]A",
            "T[T>C]C",
            "T[T>C]G",
            "T[T>C]T",
            "A[T>G]A",
            "A[T>G]C",
            "A[T>G]G",
            "A[T>G]T",
            "C[T>G]A",
            "C[T>G]C",
            "C[T>G]G",
            "C[T>G]T",
            "G[T>G]A",
            "G[T>G]C",
            "G[T>G]G",
            "G[T>G]T",
            "T[T>G]A",
            "T[T>G]C",
            "T[T>G]G",
            "T[T>G]T",
        ]
        if revcc:
            seq96 = [rev_comp(x) for x in seq96]
            print(seq96)
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [1 / 256, 1 / 256, 1 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [203 / 256, 202 / 256, 202 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [236 / 256, 199 / 256, 197 / 256],
        ]
        xlabels = [seq[0] + seq[2] + seq[6] for seq in seq96]
        i = 0

        x = 0.043
        y3 = 0.87
        y = int(ymax * 1.25)
        y2 = y + 2
        for i in range(0, 6, 1):
            panel1.add_patch(
                plt.Rectangle(
                    (x, y3),
                    0.15,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            x += 0.159

        yText = y3 + 0.06
        plt.text(
            0.1,
            yText,
            "C>A" if not revcc else "G>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.255,
            yText,
            "C>G" if not revcc else "G>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.415,
            yText,
            "C>T" if not revcc else "G>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.575,
            yText,
            "T>A" if not revcc else "A>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.735,
            yText,
            "T>C" if not revcc else "A>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.89,
            yText,
            "T>G" if not revcc else "A>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )

        if y <= 4:
            y += 4

        while y % 4 != 0:
            y += 1
        # ytick_offest = int(y/4)
        y = ymax / 1.025
        ytick_offest = float(y / 3)
        labs = np.arange(0.375, 96.375, 1)

        panel1.set_xlim([0, 96])
        # panel1.set_ylim([0, y])
        panel1.set_xticks(labs)
        # panel1.set_yticks(ylabs)
        count = 0
        m = 0
        for i in range(0, 96, 1):
            plt.text(
                i / 101 + 0.0415,
                0.02,
                xlabels[i][0],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.044,
                xlabels[i][1],
                fontsize=30,
                color=colors[m],
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                fontweight="bold",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.071,
                xlabels[i][2],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            count += 1
            if count == 16:
                count = 0
                m += 1

        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")

        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=False,
            left=True,
            labelleft=True,
            right=True,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="lightgray",
            width=2,
        )

        [i.set_color("black") for i in plt.gca().get_yticklabels()]
        if return_plot_template == False:
            pickle.dump(plot1, open(path, "wb"))
        else:
            pickle.dump(plot1, open(path, "wb"))
            return plot1


def plotSBS96(
    matrix_path,
    output_path,
    project,
    percentage=False,
    custom_text_upper=None,
    custom_text_middle=None,
    custom_text_bottom=None,
    savefig_format="pdf",
    volume=None,
    dpi=100,
    revcc=False,
):
    """Use an input matrix to create a SBS plot.

    Args:
            matrix_path: The path to a text file or a pandas DataFrame.
            output_path: Path to a directory for saving the output.
            project: Name of unique sample set
            savefig_format: Format of the output plot (pdf, png, or PIL_Image)
            volume: Path to the .pkl file containing the plot template. For Docker.
    Returns:
            Plot of the given input matrix.
    """
    plot_custom_text = False
    sig_probs = False
    pcawg = False
    plot_type = '96'

    # load custom fonts for plotting
    load_custom_fonts()

    # create the output directory if it doesn't exist
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    try:
        data = process_input(matrix_path, plot_type)
        data = reindex_sbs96(data)
        # print(data)
        sample_count = 0

        buf = io.BytesIO()
        fig_orig = make_pickle_file_SBS96(
            context="SBS96", return_plot_template=True, volume=volume, revcc=revcc,
        )
        pickle.dump(fig_orig, buf)

        figs = {}
        buff_list = {}
        ctx = data.index  # [seq[0]+seq[2]+seq[6] for seq in data.index]
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [1 / 256, 1 / 256, 1 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [203 / 256, 202 / 256, 202 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [236 / 256, 199 / 256, 197 / 256],
        ]
        colorsall = [
            [colors[j] for i in range(int(len(ctx) / 6))] for j in range(6)
        ]
        colors_flat_list = [item for sublist in colorsall for item in sublist]

        for sample in data.columns:
            buf.seek(0)
            figs[sample] = pickle.load(buf)
            panel1 = figs[sample].axes[0]

            # total_count = np.sum(data[sample].values)
            total_count = 1
            x = 0.4
            ymax = 0
            i = 0
            muts = data[sample].values
            
            if percentage:
                if total_count > 0:
                    plt.bar(
                        np.arange(len(ctx)) + x,
                        muts / total_count * 100,
                        width=0.4,
                        color=colors_flat_list,
                        align="center",
                        zorder=1000,
                    )
                    ymax = np.max(muts / total_count * 100)
                    print(ymax)
                sig_probs = True
            else:
                plt.bar(
                    np.arange(len(ctx)) + x,
                    muts,
                    width=0.4,
                    color=colors_flat_list,
                    align="center",
                    zorder=1000,
                )
                ymax = np.max(muts)
            
            x = 0.043
            y3 = 0.87
            y = int(ymax * 1.25)
            y2 = y + 2

            if y <= 4:
                y += 4

            while y % 4 != 0:
                y += 1

            y = ymax / 1.025
            ytick_offest = float(y / 3)

            # CUSTOM
            ytick_offest = 115

            if percentage:
                ylabs = [
                    0,
                    round(ytick_offest, 1),
                    round(ytick_offest * 2, 1),
                    round(ytick_offest * 3, 1),
                    round(ytick_offest * 4, 1),
                ]
                ylabels = [
                    str(0),
                    str(round(ytick_offest/100, 2)) + "%",
                    str(round(ytick_offest/100 * 2, 2)) + "%",
                    str(round(ytick_offest/100 * 3, 2)) + "%",
                    str(round(ytick_offest/100 * 4, 2)) + "%",
                ]
            else:
                ylabs = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    8,
                ]
                ylabels = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
                ]

            labs = np.arange(0.375, 96.375, 1)

            font_label_size = 30
            if not percentage:
                if int(ylabels[3]) >= 1000:
                    font_label_size = 20

            if percentage:
                if len(ylabels) > 2:
                    font_label_size = 20

            if not percentage:
                ylabels = getylabels(ylabels)

            panel1.set_xlim([0, 96])
            panel1.set_ylim([0, y])
            panel1.set_yticks(ylabs)
            if sig_probs:
                plt.text(
                    0.045,
                    0.75,
                    sample,
                    fontsize=30,
                    weight="bold",
                    color="black",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    0.045,
                    0.75,
                    sample + ": " + "{:,}".format(int(total_count)) + " subs",
                    fontsize=60,
                    weight="bold",
                    color="black",
                    fontname="Arial",
                    transform=plt.gcf().transFigure,
                )

            panel1.set_yticklabels(ylabels, fontsize=font_label_size)
            plt.gca().yaxis.grid(True)
            plt.gca().grid(
                which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1
            )
            panel1.set_xlabel("")
            panel1.set_ylabel("")

            custom_text_upper_plot = ""
            try:
                custom_text_upper[sample_count]
            except:
                custom_text_upper = False
            try:
                custom_text_middle[sample_count]
            except:
                custom_text_middle = False
            try:
                custom_text_bottom[sample_count]
            except:
                custom_text_bottom = False

            if custom_text_upper:
                plot_custom_text = True
                if len(custom_text_upper[sample_count]) > 40:
                    print(
                        "To add a custom text, please limit the string to <40 characters including spaces."
                    )
                    plot_custom_text = False
            if custom_text_middle:
                if len(custom_text_middle[sample_count]) > 40:
                    print(
                        "To add a custom text, please limit the string to <40 characters including spaces."
                    )
                    plot_custom_text = False

            if plot_custom_text:
                x_pos_custom = 0.98
                if custom_text_upper and custom_text_middle:
                    custom_text_upper_plot = (
                        custom_text_upper[sample_count]
                        + "\n"
                        + custom_text_middle[sample_count]
                    )
                    if custom_text_bottom:
                        custom_text_upper_plot += (
                            "\n" + custom_text_bottom[sample_count]
                        )

                if custom_text_upper and not custom_text_middle:
                    custom_text_upper_plot = custom_text_upper[sample_count]
                    panel1.text(
                        x_pos_custom,
                        0.78,
                        custom_text_upper_plot,
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                        ha="right",
                    )

                elif custom_text_upper and custom_text_middle:
                    if not custom_text_bottom:
                        panel1.text(
                            x_pos_custom,
                            0.72,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )
                    else:
                        panel1.text(
                            x_pos_custom,
                            0.68,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                elif not custom_text_upper and custom_text_middle:
                    custom_text_upper_plot = custom_text_middle[sample_count]
                    panel1.text(
                        x_pos_custom,
                        0.78,
                        custom_text_upper_plot,
                        fontsize=40,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                        ha="right",
                    )

            if percentage:
                plt.ylabel(
                    "Percentage of Single Base Substitutions",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                )
            else:
                plt.ylabel(
                    "Number of Single Base Substitutions",
                    fontsize=35,
                    fontname="Times New Roman",
                    weight="bold",
                )

            panel1.tick_params(
                axis="both",
                which="both",
                bottom=False,
                labelbottom=False,
                left=True,
                labelleft=True,
                right=True,
                labelright=False,
                top=False,
                labeltop=False,
                direction="in",
                length=25,
                colors="lightgray",
                width=2,
            )

            [i.set_color("black") for i in plt.gca().get_yticklabels()]
            sample_count += 1

        return output_results(
            savefig_format, output_path, project, figs, "SBS_96", dpi=dpi
        )
    except:
        print("There may be an issue with the formatting of your matrix file.")
        pdf_path = output_path + "SBS_96_plots_" + project + ".pdf"
        if os.path.isfile(pdf_path):
            os.remove(pdf_path)




spectra = pd.read_csv("../data/new_dataset/MutSpecVertebrates192.csv.gz")
spectra = spectra[spectra.Gene == 'Cytb'] # one gene
spectra_mean = spectra.groupby('Mut').MutSpec.mean().reset_index()
spectra_mean['Mut'] = spectra_mean.Mut.apply(rev_comp)

# d = pd.read_csv("../data/MutSpecALLvert.csv")
d = spectra_mean
sample_col = 'Vertebrates average spectrum'
d.columns = ['MutationType', sample_col]
d[sample_col] *= 100
# d[sample_col] = (d[sample_col] * 10000).astype(int)
first_muts = ['C>A', 'C>G', 'C>T', 'T>A', 'A>G', 'T>G',]

d1 = d[d.MutationType.str.slice(2, 5).isin(first_muts)]
d1.loc[d1.MutationType.str.slice(2,5) == 'A>G', 'MutationType'] = \
    d1[d1.MutationType.str.slice(2,5) == 'A>G'].MutationType.apply(rev_comp)

# ag = d1[d1.MutationType.str.slice(2, 5) == 'A>G']
# ag['MutationType'] = ag.MutationType.apply(rev_comp)

# d1 = d1[d1.MutationType.str.slice(2, 5) != 'A>G']
# d1 = pd.concat([d1, ag])

# print(d1[d1.MutationType.str.slice(2, 5) == 'A>G'])


d2 = d[~d.MutationType.str.slice(2, 5).isin(first_muts)]
d2.loc[d2.MutationType.str.slice(2,5) == 'T>C', 'MutationType'] = \
    d2[d2.MutationType.str.slice(2,5) == 'T>C'].MutationType.apply(rev_comp)
d2['MutationType'] = d2.MutationType.apply(rev_comp)

# print(d1.sum(), d2.sum())

d1.to_csv('/tmp/test_vert1_96.txt', sep='\t', index=False)
d2.to_csv('/tmp/test_vert2_96.txt', sep='\t', index=False)

matrix_path = '/tmp/test_vert1_96.txt'
output_path = './'
project = 'vert_1'
plotSBS96(matrix_path, output_path, project, percentage=True, revcc=False)


matrix_path = '/tmp/test_vert2_96.txt'
output_path = './'
project = 'vert_2'
plotSBS96(matrix_path, output_path, project, percentage=True, revcc=True)

