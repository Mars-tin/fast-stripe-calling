import pandas as pd

from utils import *


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
    df_tad = pd.read_csv("data/tad/{}_{}kb.csv".format(cell_type, 25))
    df_tad = df_tad.loc[df_tad['type'] == "domain"]
    df_tad = df_tad.loc[df_tad['chrom'] == chrome]
    list_tad = np.asarray([df_tad['start']//resolution, df_tad['end']//resolution]).T

    results = []
    print(' Step 4: Statistical Tests...')
    for elm in new_positions:
        [st, ed, head, tail, score] = elm
        p = stat_test(mat, orientation, st, ed, head, tail, stripe_width, window_size, list_tad)
        if p < 0.2:
            results.append([st, ed, head, tail, score, 1-p/2])
    return results
