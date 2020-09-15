import os
import numpy as np
from scipy.stats import kruskal
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


def show_heat_plot(matrix, ran=(22650, 22950), vmin=0, vmax=10):
    matrix = matrix[ran[0]:ran[1]]
    sns.set()
    sns.heatmap(matrix, vmin=vmin, vmax=vmax, cmap="Reds")
    plt.savefig("line_matrix.pdf")
    plt.show()


def pick_max_positions(mat, interval=500000, distance_range=(0, 300000),
                       resolution=25000, line_width=1, window_size=5):
    """
    :param mat: (2D ndarray) line matrix
    :param resolution: (int) resolution
    :param distance_range: (int) only consider this range off the diagonal
    :param interval: (int) minimum interval between two stripe anchors
    :param line_width: (int) stripe width (# of bins)
    :param window_size: (int) size of the window for calculating enrichment
    """

    assert interval % resolution == 0
    assert distance_range[0] % resolution == 0
    assert distance_range[1] % resolution == 0

    start, end = distance_range[0] // resolution, distance_range[1] // resolution
    size = interval // resolution
    length = mat.shape[0]
    stats = np.sum(mat[:, start:end], axis=1)
    all_pos = []

    for i in range(0, length, size):
        region = stats[i: min(i + size, length)]
        idx = int(np.argmax(region) + i)

        if idx < window_size or idx >= mat.shape[0] - window_size:
            continue

        previous = stats[max(0, idx - size): idx-1]
        later = stats[idx + 2: min(idx + size + 1, length)]

        if stats[idx] > np.max(previous) and stats[idx] > np.max(later):
            check = line_neighbor_score(mat, idx, line_width=line_width,
                                        distance_range=(start, end),
                                        window_size=window_size,
                                        neighbor_trick='mean', line_trick='min')
            if np.sum(check) > 0:
                all_pos.append(idx)
    return all_pos


def line_neighbor_score(mat, idx, line_width=1, distance_range=(20, 40), window_size=5,
                        line_trick='min', neighbor_trick='mean', metric='diff'):
    """
    :math: enrichment-score = line-trick(line) - neighbor-trick(neighbor)

    :param idx: (int) stripe location index
    :param mat: (2D ndarray) line matrix
    :param distance_range: (int) only consider this range off the diagonal
    :param line_width: (int) stripe width (# of bins)
    :param window_size: (int) size of the window for calculating enrichment
    :param line_trick: (str) 'min' or 'med'
    :param neighbor_trick: (str) 'mean' or 'med'
    :param metric:  (str) 'diff' or 'ratio'
    """

    half = int(line_width // 2)
    x1, x2 = idx - half, idx - half + line_width

    new_mat = np.zeros((distance_range[1] - distance_range[0],))
    for j in range(distance_range[0], distance_range[1]):

        if j < window_size + half or j >= mat.shape[1] - window_size - half:
            continue
        y = j - distance_range[0]

        if line_trick == 'min':
            line_score = min(np.mean(mat[x1:x2, j - window_size - half:j - half]),
                             np.mean(mat[x1:x2, j + 1 + half:j + window_size + half + 1]))
        else:
            line_score = np.median(np.concatenate(
                [mat[x1:x2, j - window_size - half:j - half],
                 mat[x1:x2, j + 1 + half:j + window_size + half + 1]]
            ))

        if neighbor_trick == 'mean':
            neighbor_score = max(np.mean(mat[idx-window_size:x1, j-window_size-half:j+window_size+half+1]),
                                 np.mean(mat[x2+1:idx+window_size+1, j-window_size-half:j+window_size+half+1]))
        else:
            neighbor_score = max(np.median(mat[idx-window_size:x1, j-window_size-half:j+window_size+half+1]),
                                 np.median(mat[x2+1:idx+window_size+1, j-window_size-half:j+window_size+half+1]))

        if metric == 'diff':
            new_mat[y] = line_score - neighbor_score
        else:
            new_mat[y] = (line_score / neighbor_score - 1) if neighbor_score != 0 else line_score
    return new_mat


def find_max_slice(arr):
    maximum, head, tail = 0, 0, 0
    score, head_temp, tail_temp = 0, 0, 0
    i = 0
    while i < len(arr):
        score = score + arr[i]
        if score <= 0:
            head_temp, tail_temp = i + 1, i + 1
            score = 0
        else:
            tail_temp = i + 1
        if score > maximum:
            head, tail, maximum = head_temp, tail_temp, score
        i += 1
    if tail > head:
        while arr[tail] < 0:
            tail -= 1
    return head, tail, maximum


def merge_positions(lst, merge_range, normalized=False, tol=2):

    def _has_overlap(pair1, pair2, tolerance):
        return pair1[0] - tolerance <= pair2[1] and pair1[1] >= pair2[0] - tolerance

    def _max_merge(neighbors):
        merged = []
        for elm in neighbors:
            matched = False
            for branch in merged:
                # [st, ed, head, tail, score] = elm
                if _has_overlap([elm[2], elm[3]], [branch[2], branch[3]], tolerance=tol):
                    branch[1] = elm[1]
                    branch[2] = min(branch[2], elm[2])
                    branch[3] = max(branch[3], elm[3])
                    branch[4] = max(branch[4], elm[4])
                    matched = True
                    break
            if not matched:
                merged.append(elm)
        if len(merged) > 1:
            merged.sort(key=lambda x: x[4], reverse=True)
        return merged

    def _normalized_merge(neighbors):
        merged = []
        for elm in neighbors:
            matched = False
            for branch in merged:
                # [st, ed, head, tail, score] = elm
                if _has_overlap([elm[2], elm[3]], [branch[2], branch[3]], tolerance=tol):
                    branch[1] = elm[1]
                    branch[2] = min(branch[2], elm[2])
                    branch[3] = max(branch[3], elm[3])
                    branch[4] += elm[4]
                    matched = True
                    break
            if not matched:
                merged.append(elm)
        for branch in merged:
            branch[4] /= (branch[1]-branch[0]+1)
        if len(merged) > 1:
            merged.sort(key=lambda x: x[4], reverse=True)
        return merged

    new_lst = []
    temp = []
    for i, (idx, head, tail, score) in enumerate(lst):
        if i != len(lst) - 1 and lst[i + 1][0] - idx < merge_range:
            temp.append([idx, idx, head, tail, score])
        else:
            if len(temp) != 0:
                temp.append([idx, idx, head, tail, score])
                if normalized:
                    new_lst.extend(_normalized_merge(temp))
                else:
                    new_lst.extend(_max_merge(temp))
                temp = []
            else:
                new_lst.append([idx, idx, head, tail, score])
    return new_lst


def stat_test(mat, orientation, st, ed, head, tail, line_width, window_size, list_tad=None):
    def _has_overlap(pair1, pair2, tolerance):
        return pair1[0] - tolerance <= pair2[1] and pair1[1] >= pair2[0] - tolerance

    half = int(line_width // 2)
    x1, x2 = st - half, ed + half + 1

    local = mat[x1:x2, head:tail].flatten()
    upper = mat[x1 - window_size:x1, head:tail].flatten()
    lower = mat[x2:x2 + window_size, head:tail].flatten()
    t1, p1 = kruskal(local, upper)
    t2, p2 = kruskal(local, lower)

    if list_tad is not None:
        tad_region = None
        for tad in list_tad:
            if _has_overlap(tad, [x1, x2], 0):
                tad_region = tad
                break
        if tad_region is not None:
            symm = []
            if orientation == 'h':
                anc = tad_region[0] + tad_region[1] - (x1 + x2) // 2 - head
                for i in range(tail - head):
                    symm.extend(mat[anc - i, head + i - half:head + i + half + 1])
            else:
                anc = tad_region[0] + tad_region[1] - (x1 + x2) // 2 + head
                for i in range(tail - head):
                    symm.extend(mat[anc + i, head + i - half:head + i + half + 1])
            symm = np.array(symm)
            t3, p3 = kruskal(local, symm)
            return np.max([p1, p2, p3])

    return max(p1, p2)
