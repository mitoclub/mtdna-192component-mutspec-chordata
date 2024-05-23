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
from pymutspec.annotation import rev_comp

from sigProfilerPlotting import (
    load_custom_fonts, process_input, reindex_sbs96, 
    make_pickle_file, getylabels, output_results,
)

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
        sample_count = 0

        buf = io.BytesIO()
        fig_orig = make_pickle_file(
            context="SBS96", return_plot_template=True, volume=volume
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

            total_count = np.sum(data[sample].values)
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
                    str(round(ytick_offest, 1)) + "%",
                    str(round(ytick_offest * 2, 1)) + "%",
                    str(round(ytick_offest * 3, 1)) + "%",
                    str(round(ytick_offest * 4, 1)) + "%",
                ]
            else:
                ylabs = [
                    0,
                    ytick_offest,
                    ytick_offest * 2,
                    ytick_offest * 3,
                    ytick_offest * 4,
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
                    fontsize=60,
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





d = pd.read_csv("./data/MutSpecALLvert.csv")
d.columns = ['MutationType', 'Vert']
# d['Vert'] = (d['Vert'] * 10000).astype(int)
first_muts = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G',]

d1 = d[d.MutationType.str.slice(2, 5).isin(first_muts)]

d2 = d[~d.MutationType.str.slice(2, 5).isin(first_muts)]
d2['MutationType'] = d2.MutationType.apply(rev_comp)

d1.to_csv('./test_vert1_96.txt', sep='\t', index=False)
d2.to_csv('./test_vert2_96.txt', sep='\t', index=False)
matrix_path = './test_vert1_96.txt'
output_path = './'
project = 'vert_1'


plotSBS96(matrix_path, output_path, project, percentage=True)
matrix_path = './test_vert2_96.txt'
output_path = './'
project = 'vert_2'
plotSBS96(matrix_path, output_path, project, percentage=True)

