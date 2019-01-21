from BITS.utils import run_command

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Aggregate all motifs under CWD and filter by some thresholds")
    parser.add_argument('motifs_fofn',
                        type=str,
                        help="File of *.motif_summary.csv file names")
    parser.add_argument('-u',
                        '--keep_modified_base',
                        action='store_true',
                        help="Not filter motifs whose methylation types are 'modified_base' [False]")
    parser.add_argument('-f',
                        '--min_frac',
                        type=float,
                        default=0.6,
                        help="Minimum methylation 'fraction'")
    parser.add_argument('-n',
                        '--min_freq',
                        type=int,
                        default=5,
                        help="Minimum 'nDetected'")
    parser.add_argument('-i',
                        '--min_ipdr',
                        type=float,
                        default=2.5,
                        help="Minimum 'meanIpdRatio'")
    parser.add_argument('-c',
                        '--min_cov',
                        type=float,
                        default=40.,
                        help="Minimum 'meanCoverage'")
    args = parser.parse_args()

    # Merge all motifs with original bin names at the first column
    run_command(f"cat {args.motifs_fofn} | while read DATA; do cat ${{DATA}} | while read line; do echo -e \"${{DATA##*/}}\t${{line}}\" | sed 's/\"//g' | sed 's/,/\t/g'; done; done | grep -v 'motifString' > all.motifs")

    # Filter the motifs using the given thresholds
    with open("all.motifs.filtered", 'w') as f:
        with open("all.motifs", 'r') as g:
            for line in g:
                line = line.strip()
                data = line.split('\t')
                if (((not args.keep_modified_base) and data[3] == "modified_base")
                    or float(data[4]) < args.min_frac
                    or int(data[5]) < args.min_freq
                    or float(data[10]) < args.min_ipdr
                    or float(data[11]) < args.min_cov):
                    continue
                f.write(f"{line}\n")
