import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

"""
Bland-Altman plot for comparing methylation-based and scRNA-based cell proportion estimates.

X-axis: mean of paired measurements (methylation + scRNA) / 2
Y-axis: difference between measurements (methylation - scRNA)
Horizontal lines: mean difference (bias) and ±1.96*SD (limits of agreement)
"""

if __name__ == '__main__':
    # --- Load and prepare data (same pipeline as plot_contrasting_bioheart_proportions.py) ---
    meth_proportions = pd.read_csv("./bioheart_methylations/Illuminas_deconv_output.csv")
    scRNA_proportions = pd.read_csv("./bioheart_methylations/CT_Longread_proportions.csv")

    meth_proportions.columns = list(meth_proportions.columns[:1]) + [c.split(".")[0] for c in meth_proportions.columns[1:]]
    meth_proportions = meth_proportions.iloc[:, :-1]

    scRNA_proportions = scRNA_proportions.groupby('FS_ID')
    pivoted = []
    for fs_id, group in scRNA_proportions:
        pivoted_df = group.pivot(index='FS_ID', columns='Cell Type', values=['Percentage']).reset_index()
        pivoted.append(pivoted_df)
    scRNA_proportions = pd.concat(pivoted, ignore_index=True)

    scRNA_proportions = scRNA_proportions.transpose()
    scRNA_proportions.columns = scRNA_proportions.iloc[0]
    scRNA_proportions = scRNA_proportions.drop(scRNA_proportions.index[0])
    scRNA_proportions.index = scRNA_proportions.index.map(lambda x: x[1])
    scRNA_proportions = scRNA_proportions.apply(lambda x: x.str.rstrip('%').astype('float') / 100)

    meth_proportions = meth_proportions.drop(meth_proportions.index[5:])
    meth_proportions.index = meth_proportions['Unnamed: 0']
    meth_proportions = meth_proportions.drop(columns=['Unnamed: 0'])

    meth_proportions_name_mapping = {
        'Monocytes_EPIC': 'Monocytes',
        'B-cells_EPIC': 'B cells',
        'NK-cells_EPIC': 'NK cells',
        'CD4T-cells_EPIC': 'CD4 T cells',
        'CD8T-cells_EPIC': 'CD8 T cells',
    }
    scRNA_proportions_name_mapping = {
        'Monos': 'Monocytes',
        'B cells': 'B cells',
        'NK cells': 'NK cells',
        'CD4 cells': 'CD4 T cells',
        'CD8 cells': 'CD8 T cells',
    }

    meth_proportions.index = meth_proportions.index.map(meth_proportions_name_mapping)
    scRNA_proportions.index = scRNA_proportions.index.map(scRNA_proportions_name_mapping)

    merged_data = pd.concat([meth_proportions.stack(), scRNA_proportions.stack()], axis=1)
    merged_data.columns = ['Methylation', 'scRNA']

    # --- Compute Bland-Altman metrics ---
    merged_data['Mean'] = (merged_data['Methylation'] + merged_data['scRNA']) / 2
    merged_data['Difference'] = merged_data['scRNA'] - merged_data['Methylation']

    mean_diff = merged_data['Difference'].mean()
    std_diff = merged_data['Difference'].std(ddof=1)
    upper_loa = mean_diff + 1.96 * std_diff
    lower_loa = mean_diff - 1.96 * std_diff

    # --- Plot ---
    sns.set_style('whitegrid')
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.set(font_scale=1.2)

    merged_data['Cell_Type'] = merged_data.index.get_level_values(0)

    cell_type_colors = {
        'Monocytes': '#1f77b4',    # blue
        'CD4 T cells': '#2ca02c',  # green
        'CD8 T cells': '#9467bd',  # purple
        'NK cells': '#d62728',     # red
        'B cells': '#ff7f0e',      # orange
    }

    scatter = sns.scatterplot(
        data=merged_data, x='Mean', y='Difference', hue='Cell_Type',
        palette=cell_type_colors, s=100, edgecolor='k', linewidth=0.5,
    )

    # Limits of agreement
    ax.axhline(y=mean_diff, color='black', linestyle='-', linewidth=1.5,
               label=f'Mean diff (bias): {mean_diff:.3f}')
    ax.axhline(y=upper_loa, color='gray', linestyle='--', linewidth=1.2,
               label=f'+1.96 SD: {upper_loa:.3f}')
    ax.axhline(y=lower_loa, color='gray', linestyle='--', linewidth=1.2,
               label=f'-1.96 SD: {lower_loa:.3f}')
    ax.axhline(y=0, color='red', linestyle=':', linewidth=1.0, alpha=0.5)

    ax.set_xlabel('Mean Proportion Estimate', fontsize=14)
    ax.set_ylabel('scRNA - Methylation', fontsize=14)
    ax.set_title('Bland-Altman Plot: Methylation vs scRNA Proportions', fontsize=16)
    ax.legend(title='Cell Type')
    ax.tick_params(axis='both', labelsize=12)

    plt.tight_layout()
    # plt.savefig('./bioheart_methylations/scrna_meth_bland_altman.png', dpi=300, bbox_inches='tight')
    plt.show()
