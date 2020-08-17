# NaiveStripe
Python script for calling stripes from HiC contact maps. This is a naive version.
- Input: `.hic` file
- Output: a `.bedpe` 2D annotation and `.signal` file

## Requirements
- `python3`
  - Packages:` numpy`, `scipy`, `matplotlib`, `seaborn`

## Usage

### Stripe Calling

```python
python PyStripe.py -i data/4943.hic -o stripes/4943
```


- `-i`, `--input`: `.hic` file path
- `-o`, `--output`: output file name

### Stripe Comparison

```python
python compare.py
```

### Average Stripe Analysis

```python
python average.py [-p]
```
