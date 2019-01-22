A bundle of the custom codes used in the paper: XXX.

It contains three modules as follows:

## 1. Filtering of human-derived reads from .bax.h5 files of PacBio raw data (`human_reads_removal/`)

[DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB), [DALIGNER](https://github.com/thegenemyers/DALIGNER), [a fork of DALIGNER](https://github.com/PacificBiosciences/DALIGNER), [DAMAPPER](https://github.com/thegenemyers/DAMAPPER) must be installed in advance.

Reads that are probably derived from the host human genome are removed. First, subreads are mapped to a human reference genome by `map_to_human_genome.sh`, and then `list_filter_subreads.py` outputs a list of the subread names which are covered by the human reference genome according to a given coverage threshld. Then `replace_bax_reads_as_meaningless.py` replaces the sequences of the specified subreads appearing in a specified .bax.h5 file into 'AAA...' so that anyone cannot obtain the original sequence information. In addition, using `filter_subreads.py`, one can remove the sequences from a fasta file.

## 2. Genome assembly and circular contig detection using modified FALCON and external binning information (`assembly/`)


## 3. Host bacteria prediction of extrachromosomal mobile genetic elements (eMGEs) using DNA methylation motifs (`methylation_matching/`)

