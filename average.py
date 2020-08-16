import pandas as pd
import argparse

from utils import *


def reshape_triangle(mat, size=(51, 100)):
    orig = np.ones((size[1]+size[0]-1, size[1])) * -1
    for i in range(size[1]):
        orig[i:i+size[0], i] = mat[:, i]
    return orig[:, ::-1].T


def reshape_square(mat_h, mat_v, size=(51, 100)):
    orig = np.ones((size[1]+size[0]-1, size[1]+size[0]-1)) * -1
    for i in range(size[0]):
        orig[i, i:i+size[1]] = mat_h[i, :]
        orig[i:i+size[1], i] = mat_v[size[0]-i-1, :]
    return orig


def pile_up(full_mat, cell_type, chrom, orient, size=(51, 100), resol=10000, penalized=False):
    mat, count = np.zeros(size), 0
    f = 'stripes/{}.signal'.format(cell_type)
    df = pd.read_csv(f, sep='\t')
    df = df.loc[df['chr'] == chrom]
    df = df.loc[df['orientation'] == orient]
    if len(df) > 1:
        anchors = np.asarray([df['start'], df['end']]).T
    else:
        anchors = np.asarray([df['start'].values, df['end'].values]).T
    anchors = np.mean(anchors//resol, axis=1, dtype=np.int)
    anchors = list(anchors)

    for pos in anchors:
        count += 1
        stripe = full_mat[pos - size[0]//2: pos + size[0]//2+1, :]
        if penalized:
            for i in range(size[1]):
                stripe[:, i] *= (1-8/(i+10))
        if pos-size[0]//2 < 0:
            padded = np.zeros(size)
            padded[size[0]-stripe.shape[0]:, :] = stripe
            mat += padded
        elif pos + size[0]//2+1 > full_mat.shape[0]:
            padded = np.zeros(size)
            padded[:stripe.shape[0], :] = stripe
            mat += padded
        else:
            mat += stripe

    return mat, count


def average(cell_type, reference_genome='hg38', chroms='all', size=(51, 100),
            resolution=10000, threshold=3, shape='square', penalized=False):

    # Load in the reference chrome sizes
    chrom_lengths = load_chrom_sizes(reference_genome)
    if isinstance(chroms, str) and chroms.lower() == 'all':
        chroms = chrom_lengths.keys()
    else:
        chroms = list(chroms)

    piled_h = np.zeros(size)
    piled_v = np.zeros(size)
    count_h = 0
    count_v = 0

    for chrom in chroms:
        print(f'Dumping .hic file - {chrom}')
        hic_file = 'data/{}.hic'.format(cell_type)
        txt_file = hic2txt(hic_file, chrom, resolution=resolution, output=None)
        KR_sum = load_KR_sum(txt_file)

        print(f'Loading contact maps - {chrom}')
        mat_h = txt2line(txt_file, 'horizontal', chrom_lengths[chrom], size[1]*resolution, resolution)
        mat_v = txt2line(txt_file, 'vertical', chrom_lengths[chrom], size[1]*resolution, resolution)

        mat_h = mat_h / KR_sum * 1000
        mat_v = mat_v / KR_sum * 1000
        mat_h = np.maximum(mat_h, threshold) - threshold
        mat_v = np.maximum(mat_v, threshold) - threshold

        print(f'Piling stripes - {chrom}')
        new_h, cnt_h = pile_up(mat_h, cell_type, chrom, 'horizontal', size, resolution, penalized)
        new_v, cnt_v = pile_up(mat_v, cell_type, chrom, 'vertical', size, resolution, penalized)
        piled_h += new_h
        piled_v += new_v
        count_h += cnt_h
        count_v += cnt_v

    piled_h /= count_h
    piled_v /= count_v

    print(f'Reshaping the matrix as a {shape}')
    if penalized:
        vmax = 3
    else:
        vmax = 6
    if shape == 'triangle':
        orig_h = reshape_triangle(piled_h, size)
        orig_v = reshape_triangle(piled_v, size)
        mats = {'horizontal': orig_h, 'vertical': orig_v}

        for orient, mat in mats.items():
            plt.figure(figsize=(10, 6))
            sns.heatmap(mat, square=False, cmap='Reds', vmax=vmax, vmin=0, xticklabels=[], yticklabels=[])
            plt.savefig('plot/averaged_{}_{}.pdf'.format(cell_type, orient))
            # plt.show()
    else:
        orig = reshape_square(piled_h, piled_v, size)
        plt.figure(figsize=(12, 10))
        sns.heatmap(orig, square=False, cmap='Reds', vmax=vmax, vmin=0, xticklabels=[], yticklabels=[])
        plt.savefig('plot/averaged_{}.pdf'.format(cell_type))
        # plt.show()
    print(f'Done')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-w', '--window_size',
        type=int,
        default=51,
        help='size of the window for stripe surroundings'
    )
    parser.add_argument(
        '-m', '--max_range',
        type=int,
        default=100,
        help='max distance off the diagonal'
    )
    parser.add_argument(
        '-r', '--resolution',
        type=int,
        default=10000,
        help='resolution'
    )
    parser.add_argument(
        '-t', '--threshold',
        type=int,
        default=3,
        help='threshold of filters'
    )
    parser.add_argument(
        '-s', '--shape',
        type=str,
        default='square',
        help='square or triangle'
    )
    parser.add_argument(
        '--penalized',
        dest='penalized',
        action='store_true'
    )
    parser.set_defaults(penalized=False)
    args, _ = parser.parse_known_args()

    cells = ['HSPC_downsample', '4943', '5832_CD34pos', '6281_DMSO', '6527', '5577', 'O3', 'MV411', 'Kasumi']
    for cell in cells:
        average(cell, size=(args.window_size, args.max_range),
                resolution=args.resolution, threshold=args.threshold,
                shape=args.shape, penalized=args.penalized
                )

