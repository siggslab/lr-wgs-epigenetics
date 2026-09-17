#!/usr/bin/python3 -u

# Copied and edited from the UXM github, for specific cell-type-colours

import pandas as pd
import numpy as np
import os.path as op
import math
import os
import matplotlib
if 'DISPLAY' not in os.environ.keys():
    matplotlib.use('Agg')
import matplotlib.pylab as plt
import matplotlib.cm
import matplotlib.colors
# from homog_mem import *
import argparse


# Plotting parameters:
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)          # figure size
COLOR_MAP = 'tab10'         # color map. See https://matplotlib.org/users/colormaps.html

class PlotDeconv:
    def __init__(self, args):
        self.args = args
        self.min_rate = .0
        self.outpath = None
        self.validate_params()
        self.df = pd.read_csv(args.csv, index_col=0).fillna(0)
        self.plot_res()

    def validate_params(self):

        # input path
        # validate_file(self.args.csv)

        # output path
        self.outpath = self.args.outpath
        if self.outpath is None:
            self.outpath = op.splitext(self.args.csv)[0] + '.pdf'

        # min_rate
        min_rate = self.args.min_rate
        if not 100.0 >= min_rate >= 0:
            print('Min rate must be in range [0,100]')
            return
        if 1 > min_rate > 0:
            print('WARNING: min_rate is considered in percentages [0,100], not as a fraction [0,1]')
        self.min_rate = min_rate / 100.

    def hide_small_tissues(self):
        """
        tissues with very small contribution are grouped to the 'other' category.
        :return: The DataFrame with the new category ('other'),
                 where the low-contribution tissues are set to 0.
        """
        if self.min_rate <= 0:
            return
        others = self.df[self.df < self.min_rate].sum()
        self.df[self.df < self.min_rate] = 0.0
        self.df = pd.concat([self.df.transpose(), others.rename('other')], axis=1)
        self.df = self.df.transpose()

    def gen_bars_colors_hatches(self, nr_tissues):
        """
        Generate combinations of colors and hatches for the tissues bars
        Every tissue will get a tuple of (color, hatch)
        the last tuple is for the 'other' category, and is always black with no hatch.
        :return: a list of tuples, with length == nr_tissues
        """
        matplotlib.rcParams['hatch.linewidth'] = 0.3
        hatches = [None, 'xxx', '...', 'O', '++'][:max(1, nr_tissues // 7)]

        nr_colors = int(math.ceil(nr_tissues / len(hatches)) + 1)

        # generate bars colors:
        cmap = matplotlib.cm.get_cmap(COLOR_MAP)
        norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_colors))
        colors = [cmap(norm(k)) for k in range(nr_colors)]

        def bar_tup_i(i):
            color_ind = i % nr_colors
            hatch_ind = int(i // math.ceil(nr_tissues / len(hatches)))
            return colors[color_ind], hatches[hatch_ind]

        colors_hatches_list = [bar_tup_i(i) for i in range(nr_tissues - 1)]
        return colors_hatches_list + [((0, 0, 0, 1), None)]

    def should_skip(self, i):
        if self.args.stubs is None:
            return False
        for s in self.args.stubs:
            if s.lower() in self.df.index[i].lower():
                return False
        return True

    def gray_out_absents(self, ch_list):
        for i in np.argwhere(self.df.values.max(axis=1) <= self.min_rate).flatten():
            ch_list[i] = ((0, 0, 0, 0), None)

    def plot_res(self):

        self.hide_small_tissues()
        nr_tissues, nr_samples = self.df.shape

        # generate bars colors and hatches:
        colors_hatches = self.gen_bars_colors_hatches(nr_tissues)
        if not self.args.full_legend:
            self.gray_out_absents(colors_hatches)

        plt.figure(figsize=FIG_SIZE)
        r = [i for i in range(nr_samples)]
        bottom = np.zeros(nr_samples)
        for i in range(nr_tissues):

            # skip ref cell types not in stubs
            if self.should_skip(i):
                continue

            colour_override = {
                ## Override for UXM
                'Blood-B': ('tab:orange', None),
                'Blood-Granul': ('tab:brown', None),
                'Blood-Mono+Macro': ('tab:blue', None),
                'Blood-NK': ('tab:red', None),
                'Blood-T': ('tab:green', None),

                ## Override for Nanomix
                'granulocyte': ('tab:brown', None),
                'T-cell': ('tab:green', None),
                'B-cell': ('tab:orange', None),
                'monocyte': ('tab:blue', None),
                'NK-cell': ('tab:red', None),

                'adipocyte': ((0.5490196078431373, 0.33725490196078434, 0.29411764705882354, 1.0), 'O'),
                # 'colon_fibroblasts': ((0.8901960784313725, 0.4666666666666667, 0.7607843137254902, 1.0), 'O'),
                # 'dermal_fibroblast': ((0.4980392156862745, 0.4980392156862745, 0.4980392156862745, 1.0), 'O'),
                # 'skeletal_muscle': ((0.7372549019607844, 0.7411764705882353, 0.13333333333333333, 1.0), 'O'),
                'neuron': ((0.7372549019607844, 0.7411764705882353, 0.13333333333333333, 1.0), '...'),
            }

            celltype_name = self.df.index[i]
            if not celltype_name in colour_override:
                c = colors_hatches[i][0]
                h = colors_hatches[i][1]
            else:
                c, h = colour_override[celltype_name]

            plt.bar(r, list(self.df.iloc[i, :]),
                    edgecolor='white',
                    width=0.85,
                    label=self.df.index[i],
                    bottom=bottom,
                    color=c,
                    hatch=h)
            bottom += np.array(self.df.iloc[i, :])

        # Custom x axis
        plt.xticks(r, [w[:NR_CHRS_XTICKS] for w in self.df.columns], rotation='vertical', fontsize=9)
        plt.xlabel("sample")
        plt.xlim(-.6, nr_samples - .4)

        # Add a legend and a title
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=2)
        plt.title('Deconvolution Results\n' + op.splitext(op.basename(self.outpath))[0])

        # adjust layout, save and show
        plt.tight_layout(rect=[0, 0, .83, 1])
        plt.savefig(self.outpath)
        print(f'dumped figure to {self.outpath}')
        if self.args.show:
            plt.show()


def main():
    args = parse_args()

    # validate parameters
    PlotDeconv(args)



def parse_args():
    a = [
        # './nanomix_nnls.deconv.csv',
        # './nanomix_llse.deconv.csv',
        './nanomix_mmse.deconv.csv',
        # './uxm_synth.deconv.csv',
        '--show',
        '--outpath',
        # './nanomix_nnls_deconv.png',
        # './nanomix_llse_deconv.png',
        './nanomix_mmse_deconv.png',
        # './uxm_synth_deconv.png',
    ]

    parser = argparse.ArgumentParser()
    parser.add_argument('csv', help='Deconvolution output csv to plot')
    parser.add_argument('--outpath', '-o', help='output. Default is the same name as CSV, but PDF suffix')
    parser.add_argument('--show', action='store_true', help='Show the figure in a pop up window')
    parser.add_argument('--min_rate', type=float, default=1.0,
            help='Cell types with smaller rates are combined to the "Other" category. Range 0-100 (percentages) [1]')
    parser.add_argument('--stubs', nargs='+',
            help='Show only reference cell types that match any of the stubs')
    parser.add_argument('--full_legend', action='store_true',
            help='For plotting: If set, legend show colors of all cell types' \
                    ' . Otherwise, show only for cell types found in at least one sample')
    return parser.parse_args(a)


if __name__ == '__main__':
    main()