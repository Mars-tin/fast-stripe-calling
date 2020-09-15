import argparse
import os
import numpy as np

from utils import load_chrom_sizes, hic2txt, load_KR_sum, txt2line
from callers import _stripe_caller, _deletion_caller


def caller1(caller_func, hic_file, output_file, reference_genome='hg38', chroms='all',
            resolution=25000, max_range=3000000,
            interval=500000, min_length=500000, closeness=1000000,
            stripe_width=1, merge=1, window_size=5,
            static_filtered=False, threshold=3):
    """
    :param caller_func: (func) caller function
    :param hic_file: (str) .hic file path
    :param output_file: (str) output filename bedpe path
    :param reference_genome: (str) reference genome
    :param chroms: (str or list) chromosomes to calculate, default: 'all'
    :param resolution: (int) resolution
    :param max_range: (int) only consider this range off the diagonal
    :param interval: (int) minimum interval between two stripe anchors (for accelerating calculation)
    :param min_length: (int) minimum length of stripes
    :param closeness: (int) maximum distance off the diagonal
    :param stripe_width: (int) stripe width (# of bins)
    :param merge: (int) merge stripes within this range (# of bins)
    :param window_size: (int) size of the window for calculating enrichment
    :param static_filtered: (bool) if matrix is filtered
    :param threshold: (float) threshold of filters
    """
    # Load in the reference chrome sizes
    chrom_lengths = load_chrom_sizes(reference_genome)
    if isinstance(chroms, str) and chroms.lower() == 'all':
        chroms = chrom_lengths.keys()
    else:
        chroms = list(chroms)

    bedpe_file = output_file + ".bedpe"
    signal_file = output_file + ".signal"

    f_bedpe = open(bedpe_file, 'w')
    f_signal = open(signal_file, 'w')
    f_bedpe.write('chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tenrichment\n')
    f_signal.write('chr\tstart\tend\thead\ttail\torientation\tenrichment\tpval\n')
    txt_file = None
    cell_type = hic_file.split('/')[1][:-4]

    for chrom in chroms:

        # Step 0: dump .hic file
        print(f'Dumping .hic file - {chrom}')
        txt_file = hic2txt(hic_file, chrom, resolution=resolution, output=txt_file)
        KR_sum = load_KR_sum(txt_file)

        # Step 1: obtain all horizontal and vertical lines
        print(f'Loading contact maps - {chrom}')
        mat_h = txt2line(txt_file, 'horizontal', chrom_lengths[chrom], max_range, resolution)
        mat_v = txt2line(txt_file, 'vertical', chrom_lengths[chrom], max_range, resolution)

        mat_h = mat_h / KR_sum * 1000
        mat_v = mat_v / KR_sum * 1000

        mat_h = np.maximum(mat_h, threshold) - threshold
        mat_v = np.maximum(mat_v, threshold) - threshold

        fil = 10 if static_filtered else 0

        print(f'Calculating stripes - {chrom} - horizontal')
        stripes_h = caller_func(mat_h, orientation='h', max_range=max_range, resolution=resolution,
                                interval=interval, min_length=min_length, closeness=closeness,
                                stripe_width=stripe_width, merge=merge, window_size=window_size,
                                chrome=chrom, cell_type=cell_type)

        if len(stripes_h) > 0:
            if not static_filtered:
                fil = np.quantile(np.array(stripes_h)[:, 4], 0.1)
            for elm in stripes_h:
                st, ed, head, tail, enr, pval = elm
                if enr > fil:
                    x1, x2 = st * resolution, (ed + 1) * resolution
                    y1, y2 = max((st + head) * resolution, x2), (ed + tail) * resolution  # avoid overlap
                    f_signal.write(f'{chrom}\t{x1}\t{x2}\t{y1}\t{y2}\thorizontal\t{enr}\t{pval}\n')
                    f_bedpe.write(f'{chrom}\t{x1}\t{x2}\t{chrom}\t{y1}\t{y2}\t0,255,0\t{enr}\n')  # green

        print(f'Calculating stripes - {chrom} - vertical')
        stripes_v = caller_func(mat_v, orientation='v', max_range=max_range, resolution=resolution,
                                interval=interval, min_length=min_length, closeness=closeness,
                                stripe_width=stripe_width, merge=merge, window_size=window_size,
                                chrome=chrom, cell_type=cell_type)

        if len(stripes_v) > 0:
            if not static_filtered:
                fil = np.quantile(np.array(stripes_v)[:, 4], 0.1)
            for elm in stripes_v:
                st, ed, head, tail, enr, pval = elm
                if enr > fil:
                    x1, x2 = max(0, (st - tail)) * resolution, (st - head) * resolution
                    y1, y2 = st * resolution, (ed + 1) * resolution
                    f_signal.write(f'{chrom}\t{y1}\t{y2}\t{x1}\t{x2}\tvertical\t{enr}\t{pval}\n')
                    f_bedpe.write(f'{chrom}\t{x1}\t{x2}\t{chrom}\t{y1}\t{y2}\t0,255,0\t{enr}\n')  # blue

    f_bedpe.close()
    f_signal.close()
    if txt_file is not None:
        os.remove(txt_file)


def PyStripe(args):
    chromosomes = ['chr' + elm for elm in args.chromosomes.split(',')] if args.chromosomes.lower() != 'all' else 'all'
    if args.feature == 'stripe':
        caller1(
            _stripe_caller,
            args.input,
            args.output,
            reference_genome=args.rg,
            chroms=chromosomes,
            resolution=args.resolution,
            max_range=args.max_distance,
            interval=300000,
            min_length=args.min_length,
            closeness=args.min_distance,
            stripe_width=args.width,
            merge=args.merge,
            window_size=args.window_size,
            static_filtered=args.filter == 'False'
        )
    elif args.feature == 'deletion':
        caller1(
            _deletion_caller,
            args.input,
            args.output,
            reference_genome=args.rg,
            chroms=chromosomes,
            resolution=args.resolution,
            max_range=args.max_distance,
            interval=300000,
            min_length=args.min_length,
            closeness=args.min_distance,
            stripe_width=args.width,
            merge=args.merge,
            window_size=args.window_size,
            static_filtered=args.filter == 'False'
        )
    elif args.feature == 'het-deletion':
        pass
    elif args.feature == 'inversion':
        pass
    else:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input',
        type=str,
        default=None,
        help='.hic file path'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='chr1_10kb_stripes_v2.txt',
        help='output bedpe path'
    )
    parser.add_argument(
        '-f', '--feature',
        type=str,
        default='stripe',
    )
    parser.add_argument(
        '-r', '--resolution',
        type=int,
        default=10000,
        help='resolution'
    )
    parser.add_argument(
        '--rg',
        type=str,
        default='hg38',
        help='reference genome'
    )
    parser.add_argument(
        '-c', '--chromosomes',
        type=str,
        default='all',
        help='chromosomes, separated by comma, e.g. 1,2,3. Can also be "all".'
    )
    parser.add_argument(
        '--max_distance',
        type=int,
        default=3000000,
        help='max distance off the diagonal to be calculated'
    )
    parser.add_argument(
        '--min_length',
        type=int,
        default=100000,
        help='minimum length of stripes'
    )
    parser.add_argument(
        '--min_distance',
        type=int,
        default=1500000,
        help='threshold for removing stripes too far away from the diagonal'
    )
    parser.add_argument(
        '--width',
        type=int,
        default=3,
        help='stripe width (# of bins)'
    )
    parser.add_argument(
        '--merge',
        type=int,
        default=3,
        help='merge stripes which are close to each other (# of bins)'
    )
    parser.add_argument(
        '--window_size',
        type=int,
        default=10,
        help='size of the window for calculating enrichment score'
    )
    parser.add_argument(
        '--filter',
        type=str,
        default='True',
        help=''
    )
    args, _ = parser.parse_known_args()
    PyStripe(args)
