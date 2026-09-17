import pandas as pd
import os

# directory = "./nanomix_deconvs/nnls"
directory = "./nanomix_deconvs/llse"
# directory = "./nanomix_deconvs/mmse"

if __name__ == '__main__':
    tables = []
    for f in os.listdir(directory):
        table = pd.read_table(os.path.join(directory, f))
        ## Naming stuff
        col_name = f.split('.')[0]
        s = col_name.split('_x_')
        if len(s) > 1:
            col_name = s[0] + ' w/out ' + s[1]

        table = table.rename(columns={'proportion':col_name, 'cell_type':'CellType'})
        table = table.set_index('CellType')
        tables.append(table)

    big_table = pd.concat(tables, axis=1)
    # re-order columns manually
    order = ['Myeloid', 'Neut', 'CD4', 'CD8', 'B', 'NK', 'blood', 'pbmc', 'pbmc w/out mono', 'pbmc w/out cd4', 'pbmc w/out cd8', 'pbmc w/out t', 'pbmc w/out b', 'pbmc w/out nk']
    big_table = big_table[order]

    # Also row-reorder
    # Mono, B, CD4, CD8, Neut, NK
    new_order = ['monocyte', 'B-cell', 'T-cell', 'granulocyte', 'NK-cell']
    for idx in big_table.index:
        if idx not in new_order:
            new_order.append(idx)
    big_table = big_table.loc[new_order]


    # big_table.to_csv("./nanomix_nnls.deconv.csv", sep=',', header=True)
    big_table.to_csv("./nanomix_llse.deconv.csv", sep=',', header=True)
    # big_table.to_csv("./nanomix_mmse.deconv.csv", sep=',', header=True)


