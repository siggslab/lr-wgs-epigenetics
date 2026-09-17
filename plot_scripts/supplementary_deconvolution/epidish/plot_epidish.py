import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == '__main__':
    data = pd.read_csv("~/Downloads/only_pures.deconv.csv")
    data = data.set_index('sample')
    # Re-order columns, so default colouring is correct
    data = data[['Mono', 'B', 'CD4T', 'NK', 'CD8T', 'Neutro', 'Eosino']]

    width = 0.5
    fig, ax = plt.subplots(figsize=(12, 12))
    bottom = np.zeros(len(data))

    for c, v in data.items():  # cell_type, values
        sample = v.keys()
        frac = np.array(list(v))

        p = ax.bar(sample, frac, label=c, bottom=bottom)
        bottom += frac

ax.set_title("Epidish Deconv on Purified Leukocytes")
ax.legend(loc="upper right")
# ax.tick_params("x", rotation=45, ha='left')
plt.setp( ax.xaxis.get_majorticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.show()
# plt.savefig("epidish_pures.png")