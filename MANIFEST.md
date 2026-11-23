# Deposit contents

Supporting data and code for the manuscript "Prophage mobilization competence
and immunity class structure the combinatorial carriage of Shiga toxin subtypes
in Escherichia coli".

| File | Size | Description |
|---|---|---|
| PDG000000004.6464.amr.metadata.tsv.gz | 61 MB | Complete NCBI Pathogen Detection release as downloaded, 596,782 records |
| analysis_table_PDG000000004.6464.tsv | 16 MB | Parsed per-record analysis fields, 144,707 rows |
| master.py | | Computes every value reported in the manuscript |
| figs.py | | Regenerates Figures 1 to 3 |
| analyze.py, analyze2.py | | Diagnostic scripts |
| results.json | | Output of master.py |
| README_data_provenance.txt | | Release, retrieval date, operational definitions |
| PROMPTS.md | | Record of large language model assistance |

The intermediate file stx_positive_records_PDG000000004.6464.tsv (129 MB) is not
deposited because it is reproduced exactly from the release file with:

    zcat PDG000000004.6464.amr.metadata.tsv.gz \
      | awk -F'\t' 'NR==1{print; next} $63 ~ /stx[12][a-z]?_operon/ {print}' \
      > stx_positive_records_PDG000000004.6464.tsv
