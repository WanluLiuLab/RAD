# RAD

Basic usage:

```python
python3 RAD.py \
    --species hg38 \
    --up_genes ./example_data/upregulated_genes.txt \
    --down_genes ./example_data/downregulated_genes.txt \
    --peak_file ./example_data/region.bed \
    -uid 1001 \
    --distance 10_000
```

Output file will be saved in the `./static/file` directory.