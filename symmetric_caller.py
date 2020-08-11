import numpy as np
import pandas as pd
from scipy.stats import kruskal


def _stripe_caller(mat, orientation='h',
                   max_range=3000000, resolution=25000,
                   interval=200000, filter_rate=0.95,
                   min_length=500000, closeness=1000000,
                   stripe_width=1, merge=1, window_size=5,
                   chrome='chr1', cell_type='hspc'):
    """
    :param mat: (2D ndarray) line matrix
    :param resolution: (int) resolution
    :param max_range: (int) only consider this range off the diagonal
    :param interval: (int) minimum interval between two stripe anchors
    :param min_length: (int) minimum length of stripes
    :param filter_rate: (float) percentage of filtered anchors
    :param closeness: (int) maximum distance off the diagonal
    :param stripe_width: (int) stripe width (# of bins)
    :param merge: (int) merge stripes within this range (# of bins)
    :param window_size: (int) size of the window for calculating enrichment
    """
    # Step 1: for different distance ranges pick the "local maximum" positions
    print(' Step 1: Finding local maximum for different contact distances...')
    positions = {}
    for dis in range(0, max_range - min_length + 1, interval // 2):
        print(f'  {dis}-{dis + interval}')
        distance_range = (dis, dis + interval)
        pos_h = pick_max_positions(mat, interval=interval,
                                   distance_range=distance_range,
                                   resolution=resolution, line_width=stripe_width,
                                   window_size=window_size)
        for p in pos_h:
            if p not in positions:
                positions[p] = []
            positions[p].append(distance_range)
    print(' A total of {} positions are located'.format(len(positions)))

    # Step 2: find the accurate range of stripe
    print(' Step 2: Finding the spanning range for each stripe...')
    all_positions = []
    lst = sorted(positions.keys())
    for i, idx in enumerate(lst):
        if idx <= window_size or idx >= mat.shape[0] - window_size:
            continue
        arr = line_neighbor_score(mat, idx, line_width=stripe_width,
                                  distance_range=(0, max_range // resolution),
                                  window_size=window_size,
                                  neighbor_trick='mean', line_trick='med',
                                  metric='ratio')
        head, tail, max_val = find_max_slice(arr)
        if max_val > 0:
            all_positions.append((idx, head, tail, max_val))
    print(" Succeed sliced {} locations".format(len(all_positions)))

    # Step 3: Merging and filtering
    print(' Step 3: Merging neighboring stripes...')
    all_positions = merge_positions(all_positions, merge_range=merge)

    new_positions = []
    for elm in all_positions:
        [_, _, head, tail, _] = elm
        if (tail - head) * resolution >= min_length \
                and head * resolution <= closeness:
            new_positions.append(elm)
        else:
            pass
    print('A total of {} positions are located after merging'.format(len(new_positions)))

    # Step 4: Statistical test on symmetric stripe region
    df_tad = pd.read_csv("data/{}_{}kb.csv".format(cell_type, 25))
    df_tad = df_tad.loc[df_tad['type'] == "domain"]
    df_tad = df_tad.loc[df_tad['chrom'] == chrome]
    list_tad = np.asarray([df_tad['start']//resolution, df_tad['end']//resolution]).T

    results = []
    print(' Step 4: Statistical Tests...')
    for elm in new_positions:
        [st, ed, head, tail, score] = elm
        p = stat_test(mat, list_tad, orientation, st, ed, head, tail, stripe_width, window_size)
        if p < 0.2:
            results.append([st, ed, head, tail, score, 1-p/2])
    return results


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


def stat_test(mat, list_tad, orientation, st, ed, head, tail, line_width, window_size):

    def _has_overlap(pair1, pair2, tolerance):
        return pair1[0] - tolerance <= pair2[1] and pair1[1] >= pair2[0] - tolerance

    half = int(line_width // 2)
    x1, x2 = st - half, ed + half + 1

    local = mat[x1:x2, head:tail].flatten()
    upper = mat[x1 - window_size:x1, head:tail].flatten()
    lower = mat[x2:x2 + window_size, head:tail].flatten()
    t1, p1 = kruskal(local, upper)
    t2, p2 = kruskal(local, lower)

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
