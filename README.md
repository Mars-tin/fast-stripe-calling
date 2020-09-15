# Fast HiC Feature Extraction

Python script for calling features from HiC contact maps. This is a naive version.

- Input: `.hic` file
- Output: a `.bedpe` 2D annotation and `.signal` file

## Requirements
- `python3`
  - Packages:` numpy`, `scipy`, `matplotlib`, `seaborn`

## Stripes

### Stripe Calling

```python
python FeatExtract.py -i data/4943.hic -o stripes/4943 -f stripe
```


- `-i`, `--input`: `.hic` file path
- `-o`, `--output`: output file name
- `-f`, `--feature`: the feature to extract

### Stripe Comparison

```python
python compare.py
```

### Average Stripe Analysis

```python
python average.py [-p]
```