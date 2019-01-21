from collections import defaultdict
from logzero import logger
from BITS.utils import run_command, sge_nize


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Call methylation motifs for each genome bin")
    parser.add_argument('bins_fofn',
                        type=str,
                        help="A file of fasta file paths of the bins")
    parser.add_argument('modifications_gff',
                        type=str,
                        help="modifications.gff of all contigs")
    parser.add_argument('-p',
                        '--smrtpipe_setup',
                        type=str,
                        default="/bio/package/pacbio/smrtanalysis/new/etc/setup.sh",
                        help="File path of setup.sh of SMRTpipe for MotifMaker execution")
    parser.add_argument('-j',
                        '--jar_path',
                        type=str,
                        default="/bio/package/pacbio/smrtanalysis/new/analysis/lib/java/motif-maker-0.2.one-jar.jar",
                        help="Path of MotifMaker executable jar file")
    parser.add_argument('-m',
                        '--mem_limit',
                        type=str,
                        default="15g",
                        help="Memory limit for each MotifMaker")
    parser.add_argument('-s',
                        '--min_score',
                        type=int,
                        default=30,
                        help="Min score for each MotifMaker")
    args = parser.parse_args()

    # Load header names in each bin fasta file
    bin_contig_names = defaultdict(set)
    with open(args.bins_fofn, 'r') as f:
        for bin_path in f:
            bin_path = bin_path.strip()
            with open(bin_path, 'r') as g:
                for line in g:
                    if line[0] == '>':
                        bin_contig_names[bin_path].add(line.strip()[1:])

    # Divide modifications.gff to each contig
    mod_data = defaultdict(list)
    with open(args.modifications_gff, 'r') as f:
        for line in f:
            if line[0] != '#':
                mod_data[line.split('\t')[0]].append(line.strip())

    # For each bin, aggregate and output modifications.gff, and run MotifMaker
    for bin_path, contig_names in bin_contig_names.items():
        with open(f"{bin_path}.gff", 'w') as f:
            for contig_name in contig_names:
                for line in mod_data[contig_name]:
                    f.write(f"{line}\n")
        with open(f"{bin_path}.qsub", 'w') as f:
            f.write(sge_nize('\n'.join([f". {args.smrtpipe_setup}",
                                        ' '.join([f"java -jar -Xmx{args.mem_limit} {args.jar_path} find",
                                                  f"--minScore {args.min_score}",
                                                  f"--gff {bin_path}.gff",
                                                  f"--fasta {bin_path}",
                                                  f"--output {bin_path}.motif_summary.csv"])]),
                             job_name=f"run_motifmaker",
                             n_core=1,
                             wait=False))
        run_command(f"qsub {bin_path}.qsub")
    logger.info("Submitted all MotifMaker jobs")
