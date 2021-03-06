A bundle of the custom codes used in the paper: XXX.

It contains three modules as described below.


## 1. Filtering of human-derived reads from .bax.h5 files of PacBio raw data (`human_reads_removal/`)

[DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB), [DALIGNER](https://github.com/thegenemyers/DALIGNER), [a fork of DALIGNER](https://github.com/PacificBiosciences/DALIGNER), [DAMAPPER](https://github.com/thegenemyers/DAMAPPER) and some python libraries used in the scripts must be installed in advance.

Reads that are probably derived from the host human genome are removed. First, subreads are mapped to a human reference genome by `map_to_human_genome.sh`, and then `list_filter_subreads.py` outputs a list of the subread names which are covered by the human reference genome according to a given coverage threshld. Then `replace_bax_reads_as_meaningless.py` replaces the sequences of the specified subreads appearing in a specified .bax.h5 file into 'AAA...' so that anyone cannot obtain the original sequence information. In addition, using `filter_subreads.py`, one can remove the sequences from a fasta file instead of a .bax.h5 file.


## 2. Genome assembly and circular contig detection using modified FALCON and external binning information (`assembly/`)

To obtain more accurate and reliable output sequences, we modified FALCON as follows: 1) output unitigs, the most conservative type of contigs, 2) cluster the unitigs into bins using MetaBAT, 3) map the binning information onto the assembly graph, and 4) for each node in the graph, if it has only 1 in-edge and 1 out-edge that have the same bin ID, then connect the two edges (i.e. unitigs). Finally we output the extended unitigs as contigs.

Slightly different scripts are used for generating initial unitigs (= conservative contigs consisting only of simple paths in a (bundled) string graph) (`modified_fc_ovlp_to_graph_unitig.py`) and for subsequent contigs (`modified_fc_ovlp_to_graph_contig.py`). More specifically, first of all the original FALCON is executed with PacBio subreads in order to generate all-vs-all read overlap information. Then `modified_fc_ovlp_to_graph_unitig.py` is performed with the output of LA4Falcon against the .las overlap files. After Quiver, Pilon, and MetaBAT are applied to the resulting unitig sequences, `modified_fc_ovlp_to_graph_contig.py` is performed with the overlap file and the list of bins generated by MetaBAT.


## 3. Host bacteria prediction of extrachromosomal mobile genetic elements (eMGEs) using DNA methylation motifs (`methylation_matching/`)

Before running the codes, one must finish PacBio assembly and methylation calling using all contigs and reads, and also binning of the contigs into each putative bacterium and selection of contigs of mobile genetic elements (MGEs) you want to associate with bacteria. Input files are as follows:

* `bins.fofn`
  * A file in which each line is a path to a multi-fasta file of contigs corresponding to a single bacterial genome
* `mges.fofn`
  * Same format as bins.fofn but corresponding to MGEs
* `modifications.gff`, `modifications.csv`
  * Output of SMRT pipe

The following shows how to use the scripts. Initially call methylation motifs for each genome bin.

```
$ python3 binwise_motif_calling.py bins.fofn modifications.gff
```

Now you have a `*.motif_summary.csv` file for each bin in `bins.fofn`. Then aggregate the motifs with some thresholds.

```
$ find /root/of/*.motif_summary.csv -name "*.motif_summary.csv" > motifs.fofn
$ python3 aggregate_and_filter_motifs.py motifs.fofn
```

Now you have `all.motifs` and `all.motifs.filtered` files. Next, for each motif, calculate IPD ratio for each genome bin and eMGE.

```
$ python3 calc_ipdr_matrix.py bins.fofn mges.fofn modifications.csv
```

You can parallelize the job using -n option. And now you have `ipdr_matrix.csv`. Lastly do matching.

```
$ python3 matching.py
```

Now you have `matchings`, which records the final result.


## License

The original code of `modified_fc_ovlp_to_graph_unitig.py` and `modified_fc_ovlp_to_graph_contig.py` was developed by Pacific Biosciences and licensed as described in `LICENSE.falcon`, and the modified codes retain the same license. The other codes developed by Yoshihiko Suzuki are licensed under MIT.