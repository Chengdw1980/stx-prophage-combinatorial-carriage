NCBI Pathogen Detection data underlying the JAM Stx prophage manuscript
=======================================================================

Release analysed : PDG000000004.6464
Organism group   : Escherichia coli / Shigella
Retrieved        : 23 November 2026
Source URL       : https://ftp.ncbi.nlm.nih.gov/pathogen/Results/
                   Escherichia_coli_Shigella/PDG000000004.6464/AMR/
                   PDG000000004.6464.amr.metadata.tsv

NOTE: NCBI retains only about ten recent releases on the FTP site, so this
release will become unavailable there within roughly two weeks. The archived
copy below is the durable record.

FILES
-----
PDG000000004.6464.amr.metadata.tsv.gz
    Complete release as downloaded. 596,782 records, 67 columns.

stx_positive_records_PDG000000004.6464.tsv
    The 144,707 records in which at least one stx operon was called by
    StxTyper (as integrated into AMRFinderPlus). All original columns.

analysis_table_PDG000000004.6464.tsv
    One row per stx-positive record, with the fields used in the analysis
    parsed out: assigned stx1 and stx2 subtypes (semicolon separated),
    count of unassigned operons, fully_resolved flag, assembly level,
    contig count, assembly length, total virulence gene count, and
    presence/absence of the ten accessory markers scored
    (eae ehxA espK nleB toxB katP espP iha saa subA).
    118,756 rows carry fully_resolved = 1.

master.py
    Computes every number reported in the manuscript and its tables,
    writing results.json.  Run as:
        python3 -I master.py stx_positive_records_PDG000000004.6464.tsv
figs.py
    Regenerates Figures 1, 2 and 3 as PDF and 300 dpi TIFF from
    results.json.  Edit the OUT path at the top before running.
analyze.py, analyze2.py
    Diagnostic scripts: the assembly-quality stratification and the
    comparison of analysis-set definitions.

DEFINITIONS USED
----------------
assigned operon    token matching stx[12][a-z]_operon
unassigned operon  token matching stx[12]_operon, with or without a
                   =PARTIAL, =PARTIAL_END_OF_CONTIG or =MISTRANSLATION status
fully resolved     at least one assigned subtype and no unassigned operon
marker present     the bare gene name is returned.  A hit carrying a
                   =PARTIAL, =PARTIAL_END_OF_CONTIG or =MISTRANSLATION
                   status is NOT counted as present.  This is the rule
                   used throughout the manuscript; counting partial hits
                   as present changes several marker percentages
                   substantially (for example stx2b subA 55.8 -> 72.3).
closed             asm_level = "Complete Genome"
lineage O157       "O157" present in scientific_name or serovar
