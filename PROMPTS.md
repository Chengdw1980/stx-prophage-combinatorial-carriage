# Record of large language model assistance

Manuscript: *Prophage mobilization competence and immunity class structure the
combinatorial carriage of Shiga toxin subtypes in Escherichia coli*
Author: Davis W. Cheng, Rocky Institute of Biotechnology

This file records how a large language model assistant was used, as required by
the journal's policy on artificial intelligence. All analyses were specified,
inspected and verified by the author, who takes responsibility for every
reported value and for the interpretation presented.

## Data retrieval

The assistant was instructed to obtain the Escherichia coli / Shigella metadata
release from the NCBI Pathogen Detection FTP site, and did so without
modification. The release analysed is PDG000000004.6464, retrieved
23 November 2026. No records were added, removed or altered on retrieval.

## Analysis

The assistant was instructed to implement, in Python, the analysis described in
the Materials and methods section:

1. Retain records whose virulence genotype field contains at least one Shiga
   toxin operon call.
2. Record an operon as assigned when a subtype designation of the form
   stx1x or stx2x is returned, and as unassigned when PARTIAL,
   PARTIAL_END_OF_CONTIG or MISTRANSLATION is returned in place of a subtype.
3. Define a genome as fully resolved when at least one subtype is assigned and
   no unassigned operon is present.
4. For each stx2 subtype, count genomes carrying it as their only stx2
   determinant and the proportion of those additionally carrying any stx1
   subtype, with Wald 95 percent confidence intervals.
5. Score the accessory markers as present only when the gene name is returned
   without a partial or truncation status.
6. For genomes carrying exactly one stx1 and one stx2 subtype, compare observed
   pair counts with counts expected under independent assortment.
7. Stratify co-carriage of stx2a and stx2c by assembly contig count within
   lineage.
8. Compare assembly lengths of closed O157 genomes carrying both stx2a and
   stx2c against those in which stx2a is the sole toxin determinant, using the
   Mann-Whitney U test with a normal approximation and tie correction.

The scripts produced are master.py (all reported values) and figs.py
(Figures 1 to 3). analyze.py and analyze2.py are diagnostic scripts.

## Verification performed at the author's direction

The assistant was instructed to test whether the difference in accessory gene
content between the combining and solitary classes could be an artefact of
assembly quality. The comparison was repeated within strata of assembly contig
count and the difference persisted in every stratum, with marker prevalence
rising rather than falling as assemblies fragmented. This result is reported in
the Results section.

The assistant was also instructed to check the internal consistency of the
manuscript against the data. This identified that the accessory-content
analysis had been computed over genomes carrying a single stx2 subtype without
the stx1 exclusion stated in an earlier draft, that "closed" denoted an
assembly level of complete genome, and that markers had been scored without
counting partial hits. The Materials and methods section was corrected to
describe the analysis as performed.

## Drafting

The assistant assisted in drafting and revising the text. The scientific
content, the direction of the work and the interpretation are the author's.

## Reproduction

    python3 -I master.py stx_positive_records_PDG000000004.6464.tsv
    python3 -I figs.py

The first writes results.json containing every reported value; the second
regenerates Figures 1 to 3.
