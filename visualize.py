import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import numpy as np


def get_size(ch):
    sizes = {
        'chr1': 248956422,
        'chr2': 242193529,
        'chr3': 198295559,
        'chr4': 190214555,
        'chr5': 181538259,
        'chr6': 170805979,
        'chr7': 159345973,
        'chr8': 145138636,
        'chr9': 138394717,
        'chr10': 133797422,
        'chr11': 135086622,
        'chr12': 133275309,
        'chr13': 114364328,
        'chr14': 107043718,
        'chr15': 101991189,
        'chr16': 90338345,
        'chr17': 83257441,
        'chr18': 80373285,
        'chr19': 58617616,
        'chr20': 64444167,
        'chr21': 46709983,
        'chr22': 50818468,
        'chrX': 156040895,
        'chrY': 57227415
    }
    if ch == 'all':
        return sum(sizes.values())
    return sizes[ch]


def has_intersection(domain1, domain2):
    return domain1[0] <= domain2[1] and domain1[1] >= domain2[0]


def comp(cell1, cell2, tolerance=10000):
    cell1_file = "stripes/{}.signal".format(cell1)
    cell2_file = "stripes/{}.signal".format(cell2)

    unique_cell2 = []
    unique_cell1 = []
    shared = []

    chromes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
               'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
               'chr21', 'chr22', 'chrX', 'chrY']

    df_cell1 = pd.read_csv(cell1_file, sep='\t')
    df_cell2 = pd.read_csv(cell2_file, sep='\t')

    for ch in chromes:
        df_cell1_ch = df_cell1.loc[df_cell1['chr'] == ch]
        df_cell2_ch = df_cell2.loc[df_cell2['chr'] == ch]

        for direction in ["horizontal", "vertical"]:
            df_cell1_dir = df_cell1_ch.loc[df_cell1_ch['orientation'] == direction]
            df_cell2_dir = df_cell2_ch.loc[df_cell2_ch['orientation'] == direction]

            if len(df_cell1_dir) > 1:
                anc_cell1 = list(np.asarray([df_cell1_dir['start'], df_cell1_dir['end']]).T)
            else:
                anc_cell1 = list(np.asarray([df_cell1_dir['start'].values, df_cell1_dir['end'].values]).T)

            if len(df_cell2_dir) > 1:
                anc_cell2 = list(np.asarray([df_cell2_dir['start'], df_cell2_dir['end']]).T)
            else:
                anc_cell2 = list(np.asarray([df_cell2_dir['start'].values, df_cell2_dir['end'].values]).T)

            if len(anc_cell2) == 0 or len(anc_cell1) == 0:
                continue

            for anchor_cell2 in anc_cell2:
                anchor_cell2[0] -= tolerance / 2
                anchor_cell2[1] += tolerance / 2
            for anchor_cell1 in anc_cell1:
                anchor_cell1[0] -= tolerance / 2
                anchor_cell1[1] += tolerance / 2

            while len(anc_cell1) > 0 and len(anc_cell2) > 0:
                s1 = anc_cell1[0]
                s2 = anc_cell2[0]
                if has_intersection(s2, s1):
                    df_domain1 = df_cell1_dir.loc[df_cell1_dir['start'] == s1[0] + tolerance/2]
                    df_domain2 = df_cell2_dir.loc[df_cell2_dir['start'] == s2[0] + tolerance/2]
                    d1 = [df_domain1['head'].values[0], df_domain1['tail'].values[0]]
                    d2 = [df_domain2['head'].values[0], df_domain2['tail'].values[0]]
                    if has_intersection(d1, d2):
                        shared.append([min(s2[0], s1[0]), max(s2[1], s1[1])])
                    else:
                        unique_cell1.append(s1)
                        unique_cell2.append(s2)
                    anc_cell2.pop(0)
                    anc_cell1.pop(0)
                elif s2[1] < s1[0]:
                    unique_cell2.append(s2)
                    anc_cell2.pop(0)
                else:
                    unique_cell1.append(s1)
                    anc_cell1.pop(0)
            for anchor_cell2 in anc_cell2:
                unique_cell2.append(anchor_cell2)
            for anchor_cell1 in anc_cell1:
                unique_cell1.append(anchor_cell1)

    print("shared: {}".format(len(shared)))
    print("{}: {}".format(cell1, len(unique_cell1)))
    print("{}: {}".format(cell2, len(unique_cell2)))
    return len(unique_cell1), len(unique_cell2), len(shared)


def venn(cell1, cell2, resolution=10, save=False):
    plt.figure()
    domain = comp(cell1, cell2, tolerance=resolution*1000)
    venn2(subsets=domain, set_labels=[cell1, cell2])
    plt.title('Stripe Comparison at {}kb'.format(resolution))
    plt.tight_layout()
    if save:
        plt.savefig('out/{}_{}.pdf'.format(cell1, cell2), format='pdf')
    else:
        plt.show()


if __name__ == "__main__":
    # cells = ['HSPC_downsample', '4943', '5832_CD34pos', '6281_DMSO', '6527', '5577', 'O3', 'MV411', 'Kasumi']
    cells = ['HSPC_downsample', '4943', '5832_CD34pos', '6527', '5577', 'O3', 'MV411', 'Kasumi']
    for i in range(len(cells)):
        for j in range(i+1, len(cells)):
            venn(cells[i], cells[j], save=True)
