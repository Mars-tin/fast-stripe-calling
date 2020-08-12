# PyStripe
Python script for calling stripes from HiC contact maps.
- Input: .hic file
- Output: .bedpe 3D annotation and .signal file

## Requirements
- python3
- Packages: numpy & scipy

## Usage

```config
 >>> python PyStripe.py -i data/4943.hic -o stripes/4943
```


- '-i', '--input': .hic file path
- '-o', '--output': output bedpe and signal name
