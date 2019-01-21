import pandas as pd
from logzero import logger
from BITS.utils import run_command


def methylated_motifs(ipdr_matrix, name, min_ipdr):
    ipdr_vector = ipdr_matrix.loc[name]
    # NOTE: NaN means the value was filtered by the threshold, and -1. means the motif was absent in the target
    return (set(ipdr_vector[ipdr_vector >= min_ipdr].index),   # methylated motifs
            set(ipdr_vector[ipdr_vector == -1.].index))   # absent motifs


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Do matching between bacterial bins and MGEs")
    parser.add_argument('bins_fofn',
                        type=str,
                        help="File of file names of bin fasta files")
    parser.add_argument('mges_fofn',
                        type=str,
                        help="File of file names of MGE fasta files")
    parser.add_argument('-m',
                        '--ipdr_matrix',
                        type=str,
                        default="ipdr_matrix.csv",
                        help="Mean IPD ratio for ecah target and for each methylation motif")
    parser.add_argument('-i',
                        '--min_ipdr',
                        type=float,
                        default=2.5,
                        help="Minimum 'meanIpdRatio'")
    args = parser.parse_args()

    # Load matrix and filter small IPD ratio values
    ipdr_matrix = pd.read_table(args.ipdr_matrix, sep='\t')

    # Load bin names and MGE names
    with open(args.bins_fofn, 'r') as f:
        bin_names = [t_path.strip().split('/')[-1] for t_path in f]
    with open(args.mges_fofn, 'r') as f:
        mge_names = [t_path.strip().split('/')[-1] for t_path in f]

    # Do matching
    matchings = pd.DataFrame()
    for mge_name in mge_names:
        mge_met, mge_na = methylated_motifs(ipdr_matrix, mge_name, args.min_ipdr)
        logger.debug(f"MGE {mge_name}:\nmet: {mge_met}\nna: {mge_na}")
        if len(mge_met) == 0:
            continue
        for bin_name in bin_names:
            bin_met, bin_na = methylated_motifs(ipdr_matrix, bin_name, args.min_ipdr)
            logger.debug(f"BIN {bin_name}:\nmet: {bin_met}\nna: {bin_na}")
            if (mge_met - bin_na) == (bin_met - mge_na):   # exact matching except NaN
                logger.debug(f"MATCH between {mge_name} and {bin_name}")
                matchings = matchings.append([[mge_name, bin_name, list(mge_met), list(bin_met)]])
    if len(matchings) > 0:
        matchings.columns = ("MGE", "bin", "MGE_motifs", "bin_motifs")
        matchings.to_csv("matchings", sep='\t')
    else:
        run_command("touch matchings")
