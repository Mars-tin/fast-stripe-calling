import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def load_chrom_sizes(reference_genome):
    """
    Load chromosome sizes for a reference genome
    """
    my_path = os.path.abspath(os.path.dirname(__file__))
    f = open(os.path.join(my_path, 'data', reference_genome + '.chrom.sizes'))
    lengths = {}
    for line in f:
        [ch, l] = line.strip().split()
        lengths[ch] = int(l)
    return lengths


def hic2txt(hic_file, ch, resolution=25000, output=None):
    """
    Dump .hic file into contact lists
    :param hic_file: (str) .hic file path
    :param ch: (str) chromosome
    :param resolution: (int) resolution to use
    :param output: (str) temporary output path
    """
    if not output:
        output = "{}_{}_resol_{}.txt".format(hic_file[:-4], ch, resolution)
        if os.path.isfile(output):
            return output
    juicer = 'juicer_tools_1.11.04_jcuda.0.8.jar'
    cmd = f'java -jar {juicer} dump observed KR {hic_file} {ch} {ch} BP {resolution} {output}'
    os.system(cmd)
    return output


def txt2matrix(txt, length, resolution=25000):
    """
    Convert contact list into a numpy matrix, size: N * N

    :param txt: str, path of input .txt file
    :param length: chromosome length
    :param resolution: int, default: 25000
    """
    f = open(txt)
    n_bins = length // resolution + 1
    mat = np.zeros((n_bins, n_bins))

    for line in f:
        p1, p2, v = line.strip().split()
        if v == 'NaN':
            continue
        p1, p2, v = int(p1), int(p2), float(v)

        if max(p1, p2) >= n_bins * resolution:
            continue

        mat[p1 // resolution, p2 // resolution] += v
        if p1 // resolution != p2 // resolution:
            mat[p2 // resolution, p1 // resolution] += v

    return mat


def txt2line(txt, direction, length, max_range, resolution=25000):
    """
    :param txt: str, path of input .txt file
    :param direction: str, 'vertical' or 'horizontal'
    :param length: chromosome length
    :param max_range: int, max distance
    :param resolution: int, default: 25000
    """
    assert max_range % resolution == 0
    f = open(txt)
    n_bins = length // resolution + 1
    rg = max_range // resolution
    mat = np.zeros((n_bins, rg))
    for line in f:
        p1, p2, v = line.strip().split()
        if v == 'NaN':
            continue
        p1, p2, v = int(p1), int(p2), float(v)
        if max(p1, p2) >= n_bins * resolution:
            continue
        if p1 > p2:
            p1, p2 = p2, p1
        p1, p2 = p1 // resolution, p2 // resolution
        if p2 - p1 >= rg:
            continue
        if direction == 'horizontal':
            mat[p1, p2 - p1] += v
        else:
            mat[p2, p2 - p1] += v
    return mat


def load_KR_sum(file):
    f = open(file)
    sums = {}
    for line in f:

        lst = line.strip().split()
        if lst[2] == 'NaN':
            continue
        p1, p2, v = int(lst[0]), int(lst[1]), float(lst[2])
        if p1 not in sums:
            sums[p1] = v
        else:
            sums[p1] += v

        if p2 != p1:
            if p2 not in sums:
                sums[p2] = v
            else:
                sums[p2] += v
    all_sums = [sums[elm] for elm in sums]
    return np.mean(all_sums)


def show_heat_plot(matrix, ran=None, vmin=0, vmax=10):
    # matrix = np.transpose(matrix)[0:100]
    # matrix = np.transpose(matrix)[0:1000]
    # matrix = matrix[11000:11300]
    matrix = matrix[22650:22950]
    sns.set()
    sns.heatmap(matrix, vmin=vmin, vmax=vmax, cmap="Reds")
    plt.savefig("line_matrix.pdf")
    plt.show()
