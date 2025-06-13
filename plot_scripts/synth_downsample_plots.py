import statistics

import matplotlib.pyplot as plt
import numpy as np

cell_types = {  # Cell types we care about. Other cell types get mapped to 'other'
    'Monocytes_EPIC': 'Monocytes',
    'B-cells_EPIC': 'B Cells',
    'NK-cells_EPIC': 'NK Cells',
    'CD4T-cells_EPIC': "CD4 T Cells",
    'CD8T-cells_EPIC': 'CD8 T Cells',
    'Neutrophils_EPIC': 'Neutrophils'
}
deconv_files = [
    './synthetic_downsample/deconv_out_1.csv',
    './synthetic_downsample/deconv_out_2.csv',
    './synthetic_downsample/deconv_out_3.csv'
]


if __name__ == '__main__':
    deconv_data = {}
    for deconv_file in deconv_files:
        f = open(deconv_file, 'r')
        header = f.readline()
        header = header.strip().split(",")[1:]
        data = {h[:-2]: {} for h in header}
        line = f.readline()
        while line != "":
            split = line.strip().split(",")
            for i, h in enumerate(header):
                data[h[:-2]][split[0]] = float(split[i+1])
            line = f.readline()
        f.close()

        new_data = {}
        for col in data:
            new_data[col] = {}
            total_other = 0
            for row in data[col]:
                if row in cell_types:
                    new_data[col][cell_types[row]] = data[col][row]
                else:
                    total_other += data[col][row]
            new_data[col]['Other Cells'] = total_other

        deconv_data[deconv_file] = new_data


    cell_types = [cell_types[c] for c in cell_types]
    cell_types.append('Other Cells')

    samples = []
    for f in deconv_data:
        for sample in deconv_data[f]:
            samples.append(sample)
        break


    data_vals = {
        c: np.zeros((3, len(samples))) # np [3=repeat-mean,max,min, 7=samples]
        for c in cell_types
    }
    for n, sample in enumerate(samples):
        # n = 0
        # sample = samples[0]
        sample_data = [deconv_data[f][sample] for f in deconv_files]
        for c in cell_types:
            # c='Monocytes'
            vals = [s[c] for s in sample_data]
            data_vals[c][0, n] = statistics.mean(vals)
            data_vals[c][1, n] = max(vals)
            data_vals[c][2, n] = min(vals)

    for cell_type in cell_types:
        # c='Monocytes'
        data = data_vals[cell_type]
        data_mean, data_max, data_min = [t.squeeze(0) for t in np.split(data, [1,2], axis=0)]

        # x = np.linspace(1, len(samples), num=len(samples))
        x = [20e6, 10e6, 5e6, 2e6, 1e6, 500e3, 200e3]
        # x = [200e3, 500e3, 1e6, 2e6, 5e6, 10e6, 20e6]
        x_labels = ["20M", "10M", "5M", "2M", "1M", "500K", "200K"]

        fig, ax = plt.subplots()
        ax.fill_between(x, data_max, data_min, alpha=0.5, linewidth=0)
        ax.plot(x, data_mean, linewidth=2)
        ax.set_xscale('log')

        ax.set(
            xlim=(100e3, 50e6),
            xticks=np.array(x),
               ylim=(0, 0.6), yticks=np.arange(start=0, stop=0.6, step=0.1)
               )
        ax.invert_xaxis()
        ax.set_xticklabels(x_labels)
        ax.set_title(f"Predicted Contribution of Downsampled {cell_type}")
        ax.set_xlabel("Total Read Count")
        ax.set_ylabel("Percentage Contribution")
        plt.savefig(f"./synthetic_downsample/{cell_type}.png")
        # plt.show()



    # plt.style.use('_mpl-gallery')

    # make data
    # np.random.seed(1)
    # x = np.linspace(0, 8, 16)
    # y1 = 3 + 4 * x / 8 + np.random.uniform(0.0, 0.5, len(x))
    # y2 = 1 + 2 * x / 8 + np.random.uniform(0.0, 0.5, len(x))
    #
    # # plot
    # fig, ax = plt.subplots()
    #
    # ax.fill_between(x, y1, y2, alpha=.5, linewidth=0)
    # ax.plot(x, (y1 + y2) / 2, linewidth=2)
    #
    # ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
    #        ylim=(0, 8), yticks=np.arange(1, 8))
    #
    # plt.show()

