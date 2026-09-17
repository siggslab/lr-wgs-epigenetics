"""
On Owen's suggestion, we are plotting:
- X axis is cell types
- Y axis is proportion

For each cell type, for a bunch of samples, the proportion of the total sample made up by the cell type.

Point is, basically, to see a cluster of B cell normals, with one that's massive, and a cluster of other cell type
normals, with one that is near the floor. For the B cell expansion leukaemia
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

cell_type_colors = {
    'Neutrophils': (140/255, 86/255, 75/255),  # #8c564b
    'Monocytes': (31/255, 119/255, 180/255),  # #1f77b4
    'CD4 T cells': (44/255, 160/255, 44/255),  # #2ca02c
    'CD8 T cells': (148/255, 103/255, 189/255),  # #9467bd
    'B cells': (255/255, 127/255, 14/255),  # #ff7f0e
    'NK cells': (214/255, 39/255, 40/255),  # #d62728
}

IDENTIFIER = "TOB0901"
# IDENTIFIER = 'm84088_240919_044341_s1.pbmeth.combin.shift.filtered'

if __name__ == '__main__':
    data = pd.read_csv("./tob_methylations/Illumina_ID_frequencies_deconv_output.csv", index_col=0)
    # data = pd.read_csv("./bioheart_methylations/Illuminas_deconv_output.csv", index_col=0)

    filt = data.drop(data.index.to_list()[6:], axis=0)

    cell_type_renames = {
        'Monocytes_EPIC': 'Monocytes',
        'B-cells_EPIC': 'B cells',
        'CD4T-cells_EPIC': 'CD4 T cells',
        'NK-cells_EPIC': 'NK cells',
        'CD8T-cells_EPIC': 'CD8 T cells',
        'Neutrophils_EPIC': 'Neutrophils',
    }
    filt = filt.rename(cell_type_renames)

    # Reorder the DataFrame based on the specified order
    ordered_cell_types = ['Neutrophils', 'Monocytes', 'CD4 T cells', 'CD8 T cells', 'B cells', 'NK cells']
    filt = filt.reindex(ordered_cell_types)

    patient_only = filt[[IDENTIFIER]]
    fig, ax = plt.subplots(figsize=(8, 10))
    p = sns.stripplot(data=patient_only.transpose(),
                  jitter=False, linewidth=1, palette=cell_type_colors,
                      size=12, marker="^")

    sns.boxplot(data=filt.transpose(), whis=1, palette=cell_type_colors)
    # sns.violinplot(data=filt.transpose())

    _ = plt.xticks(rotation=45, ha='right')


    # for i, row in enumerate(filt.index.to_list()):
    #     p.text(x=i+0.1, y=filt[IDENTIFIER][row]+0.01, s=IDENTIFIER)


    plt.show()