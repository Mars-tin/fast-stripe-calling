# PyStripe
Python script for calling stripes from HiC contact maps.
- Input: .hic file
- Output: .bedpe 3D annotation and .signal file

## Requirements
- Python 3
- Packages: numpy & scipy

## Usage

```config
 >>> python PyStripe.py -i 4943.hic -o 4943_stripes
```


- '-i', '--input': .hic file path
- '-o', '--output': output bedpe and signal name
- '-r', '--resolution': default=25000, resolution
- '--rg': default='hg38', reference genome
- '--chromosomes': default='all', chromosomes, separated by comma, e.g. 1,2,3. Can also be "all".
- '--max_distance': default=4000000, max distance off the diagonal to be calculated
- '--min_length': default=500000, minimum length of stripes
- '--min_distance': default=2000000, threshold for removing stripes too far away from the diagonal
- '--width': default=3, stripe width (# of bins)
- '--merge': default=3, merge stripes which are close to each other (# of bins)
- '--window_size': default=8, size of the window for calculating enrichment score (# of bins)
- '--threshold': default=1e-3, threshold of p values'
