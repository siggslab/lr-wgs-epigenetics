import math

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.offsetbox import TextArea
from scipy.stats import zscore, linregress
import math

"""
    meth_proportions = pd.read_csv("./bioheart_methylations/Illuminas_deconv_output.csv"):
    -----------------------
    ,FS43705317.pbmeth.combined.shifted.filtered,FS43705438.pbmeth.combined.shifted.filtered,FS43705486.pbmeth.combined.shifted.filtered,FS43705515.pbmeth.combined.shifted.filtered,FS43705553.pbmeth.combined.shifted.filtered,FS43705621.pbmeth.combined.shifted.filtered,FS43705626.pbmeth.combined.shifted.filtered,FS43705637.pbmeth.combined.shifted.filtered,FS43705655.pbmeth.combined.shifted.filtered,FS43705671.pbmeth.combined.shifted.filtered,FS43705695.pbmeth.combined.shifted.filtered,FS43705696.pbmeth.combined.shifted.filtered,FS43706193.pbmeth.combined.shifted.filtered,FS43706207.pbmeth.combined.shifted.filtered,FS43706275.pbmeth.combined.shifted.filtered,FS43706318.pbmeth.combined.shifted.filtered,FS43706327.pbmeth.combined.shifted.filtered,FS43706369.pbmeth.combined.shifted.filtered,FS43706371.pbmeth.combined.shifted.filtered,FS43706390.pbmeth.combin.shift.filtered,FS43706397.pbmeth.combin.shift.filtered,FS43706406.pbmeth.combin.shift.filtered,FS43706421.pbmeth.combin.shift.filtered,FS43706450.pbmeth.combin.shift.filtered,FS43706465.pbmeth.combin.shift.filtered,m84088_240919_044341_s1.pbmeth.combin.shift.filtered
    Monocytes_EPIC,0.170,0.253,0.161,0.013,0.213,0.234,0.151,0.349,0.136,0.139,0.188,0.288,0.143,0.151,0.152,0.196,0.166,0.137,0.090,0.258,0.068,0.172,0.209,0.392,0.153,0.152
    B-cells_EPIC,0.071,0.134,0.085,0.180,0.016,0.043,0.075,0.078,0.066,0.115,0.039,0.083,0.087,0.082,0.167,0.037,0.165,0.076,0.040,0.075,0.104,0.003,0.023,0.025,0.123,0.000

    
    scRNA_proportions = pd.read_csv("./bioheart_methylations/CT_Longread_proportions.csv")
    ---------------------------------------
    CT_ID,FS_ID,Cell Type,Count,Percentage
    CT_1856,FS43705317,NK cells,876,43.71%
    CT_1856,FS43705317,CD4 cells,449,22.41%
"""
if __name__ == '__main__':
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
    max_value = max(merged_data['Methylation'].max(), merged_data['scRNA'].max())

    # Calculate Z-scores for identifying outliers
    # merged_data['Z-Score'] = merged_data.groupby(level=0).apply(lambda x: np.sqrt(zscore(x['Methylation'])**2 + zscore(x['scRNA'])**2)).droplevel(0)
    merged_data['Distance'] = np.abs(merged_data['Methylation'] - merged_data['scRNA']) / np.sqrt(2)

    # Calculate the perpendicular distance from each point to the line y=x
    merged_data['Distance'] = np.abs(merged_data['Methylation'] - merged_data['scRNA']) / np.sqrt(2)

    # Calculate R² value using the X=Y line
    y_means = np.mean(merged_data['scRNA'])
    mean_squared_error = np.sum((merged_data['scRNA'] - y_means) ** 2)
    regression_error = np.sum((merged_data['Methylation'] - y_means) ** 2)
    r_squared = regression_error / mean_squared_error
    #
    # # Calculate R2 from a new regression line (unplotted)
    slope, intercept, r_value, p_value, std_err = linregress(merged_data['Methylation'], merged_data['scRNA'])

    angle_between_lines = math.atan((1 - slope) / (slope+1))  # In radians
    angle_between_lines = math.degrees(angle_between_lines)

    # Plot the data using a scatter plot for each cell type
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.set(font_scale=1.3)
    ax.tick_params(axis='both', labelsize=16)
    p = sns.scatterplot(data=merged_data, x='Methylation', y='scRNA', hue=merged_data.index.get_level_values(0))

    # X=Y Line
    plt.plot([0,max_value], [0, max_value], color='red', linestyle='--', label='X=Y')
    # Regression Line
    plt.plot([0, max_value], [intercept, (slope*max_value) + intercept], color='blue', linestyle='--', label='Regression')
    # Angle between the lines in the legend
    # plt.plot([], [], alpha=0, label=f"Angle: {angle_between_lines:.2f}" + u'\N{DEGREE SIGN}')

    # Label the outliers
    # for i, row in merged_data.iterrows():
    #     # if row['Z-Score'] > 3:
    #     if row['Distance'] > 0.15:
    #         plt.text(row['Methylation'], row['scRNA'], f"{row.name[1]}_{row.name[0]}", fontsize=9, color='red')

    # # Add R² value to the plot
    # plt.text(max_value - 0.15, max_value - 0.05, f'R² = {r_squared:.2f}', fontsize=12, color='blue')

    # plt.title('Contrasting Methylation and scRNA Proportions')
    plt.xlabel('LRS Derived Proportion Estimates', size=20)
    plt.ylabel('scRNA Derived Proportion Estimates', size=20)
    plt.legend(title='Cell Type')
    plt.grid(True)
    plt.show()
