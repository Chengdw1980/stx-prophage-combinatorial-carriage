# Combinatorial carriage of Shiga toxin subtypes in *Escherichia coli*

Analysis code for the manuscript:

> Cheng, D. W. *Prophage mobilization competence and immunity class structure the
> combinatorial carriage of Shiga toxin subtypes in Escherichia coli.*

**Data and archived release:** https://doi.org/10.5281/zenodo.24484329

The data are not held here. The Zenodo record above contains the complete NCBI
Pathogen Detection release analysed (PDG000000004.6464, retrieved 23 November
2026), the parsed analysis table, and these same scripts, and is the citable
archive for the manuscript.

## Scripts

| File | Purpose |
|---|---|
| `master.py` | Computes every value reported in the manuscript and its tables; writes `results.json` |
| `figs.py` | Regenerates Figures 1-3 as PDF and 300 dpi TIFF from `results.json` |
| `analyze.py` | Diagnostic: accessory gene content stratified by assembly quality |
| `analyze2.py` | Diagnostic: comparison of candidate analysis-set definitions |

## Reproducing the results

Download the release file from the Zenodo record, then:

```bash
# extract the stx-positive subset
zcat PDG000000004.6464.amr.metadata.tsv.gz \
  | awk -F'\t' 'NR==1{print; next} $63 ~ /stx[12][a-z]?_operon/ {print}' \
  > stx_positive_records_PDG000000004.6464.tsv

python3 -I master.py stx_positive_records_PDG000000004.6464.tsv
python3 -I figs.py
```

Requires Python 3 (standard library only for `master.py`); `figs.py` additionally
requires NumPy and Matplotlib. The full run takes a few seconds.

## Operational definitions

- **assigned operon** — token matching `stx[12][a-z]_operon`
- **unassigned operon** — token matching `stx[12]_operon`, with or without a
  `=PARTIAL`, `=PARTIAL_END_OF_CONTIG` or `=MISTRANSLATION` status
- **fully resolved genome** — at least one assigned subtype and no unassigned operon
- **marker present** — the bare gene name is returned; a hit carrying a partial or
  truncation status is *not* counted as present
- **closed genome** — assembly level `Complete Genome`
- **O157 lineage** — `O157` present in `scientific_name` or `serovar`

`README_data_provenance.txt` documents the release and retrieval; `PROMPTS.md`
records the large language model assistance used, as required by the journal.
