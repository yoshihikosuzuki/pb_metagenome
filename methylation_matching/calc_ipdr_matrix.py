import copy
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from logzero import logger
from BITS.seq import load_fasta

RC_MAP = dict( zip("ACGTRYWSKMDHVBNacgtrywskmdhvbn", "TGCAYRWSMKHDBVNtgcayrwsmkhdbvn") )
canonical_bases = set(['A', 'C', 'G', 'T'])
ambiguous_bases = {'N': set(canonical_bases),
                   'R': set(['A', 'G']),
                   'Y': set(['C', 'T']),
                   'W': set(['A', 'T']),
                   'S': set(['C', 'G']),
                   'K': set(['G', 'T']),
                   'M': set(['A', 'C']),
                   'D': set(['A', 'G', 'T']),
                   'H': set(['A', 'C', 'T']),
                   'V': set(['A', 'C', 'G']),
                   'B': set(['C', 'G', 'T'])}


def is_match(query, pattern):
    """
    Examine match between two strings which can have ambiguous bases.
    """
    for i in range(len(query)):
        p = pattern[i].upper()
        q = query[i].upper()
        if p in canonical_bases:
            if not q == p:
                return False
        elif p in ambiguous_bases:
            if not q in ambiguous_bases[p]:
                return False
        else:
            logger.error(f"Unrecognized base: {p}")
    return True   # match at all bases


def check_motif(mod_data, cname, cseq, mseq, pos):
    ipd_sum = 0.0
    motif_count = 0
    for i in range(0, len(cseq) - len(mseq) + 1):
        if is_match(cseq[i:i+len(mseq)], mseq):
            ipd =  mod_data[cname][i + pos]
            if ipd == -1:   # too small coverage
                continue
            else:
                ipd_sum += ipd                
                motif_count += 1
    return (ipd_sum, motif_count)


def calc_ipdr(c_data):
    bin_name, contig_names, contig_seqs, mod_f, mod_rc, motifs = c_data
    rc_seqs = ["".join([RC_MAP[c] for c in seq[::-1]]) for seq in contig_seqs]
    ipd_vector = [0.0] * len(motifs)
    counter = 0
    for data in motifs:   # for each motif
        motif_seq, pos, met_type = data
        ipd_sum = 0.0
        motif_count = 0
        for i in range(len(contig_names)):
            contig_name = contig_names[i]
            contig_seq = contig_seqs[i]
            rc_seq = rc_seqs[i]
            ipd_sum_f, motif_count_f = check_motif(mod_f, contig_name, contig_seq, motif_seq, pos)   # plus方向に走査
            ipd_sum_rc, motif_count_rc = check_motif(mod_rc, contig_name, rc_seq, motif_seq, pos)   # minus方向に走査
            ipd_sum += ipd_sum_f + ipd_sum_rc
            motif_count += motif_count_f + motif_count_rc
        if motif_count < 1:
            ipd_vector[counter] = -1
        else:
            ipd_vector[counter] = ipd_sum / motif_count
        counter += 1
    return (bin_name, ipd_vector)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate mean IPD ratio of each motif sequence pattern for each both genome bin and MGE")
    parser.add_argument('bins_fofn',
                        type=str,
                        help="File of file names of bin fasta files")
    parser.add_argument('mges_fofn',
                        type=str,
                        help="File of file names of MGE fasta files")
    parser.add_argument('modifications_csv',
                        type=str,
                        help="modifications.csv of all contigs")
    parser.add_argument('-m',
                        '--filtered_motifs',
                        type=str,
                        default="all.motifs.filtered",
                        help="Aggregated motif_summary.csv")
    parser.add_argument('-n',
                        '--n_core',
                        type=int,
                        default=1,
                        help="Degree of parallelization of the jobs")
    args = parser.parse_args()

    # Load all motifs
    motifs = set()
    with open(args.filtered_motifs, 'r') as f:
        for line in f:
            data = line.strip().split('\t')
            motifs.add((data[1], int(data[2]), data[3]))   # sequence, (0-indexed) methylated position, type
    motifs = list(motifs)

    # Load all bins and contigs in the bins
    targets = defaultdict(dict)   # {bin_name or mge_name: {contig_name: contig_seq, ...}, ...}
    with open(args.bins_fofn, 'r') as f:
        for bin_path in f:
            bin_path = bin_path.strip()
            bin_contigs = load_fasta(bin_path)
            targets[bin_path] = {c_name: c_seq for c_name, c_seq in bin_contigs.items()}
    with open(args.mges_fofn, 'r') as f:
        for mge_path in f:
            mge_path = mge_path.strip()
            mge_contigs = load_fasta(mge_path)
            targets[mge_path] = {c_name: c_seq for c_name, c_seq in mge_contigs.items()}

    # Load IPD ratio values on all positions of the contigs in <targets>
    mod_f = {c_name: np.full(len(c_seq), -1.)
             for t_name, contigs in targets.items()
             for c_name, c_seq in contigs.items()}
    mod_rc = copy.deepcopy(mod_f)
    with open(args.modifications_csv, 'r') as f:
        prev_cname = ""
        for line in f:
            if line[0] == 'r':   # header
                continue
            data = line.strip().split(',')
            contig_name = data[0].lstrip("\"").rstrip("\"")
            if contig_name not in mod_f:
                continue
            pos = int(data[1])   # start from 1
            strand = data[2]
            ipd_ratio = float(data[8])
            if strand == "0":   # plus
                mod_f[contig_name][pos - 1] = ipd_ratio
            else:   # minus
                mod_rc[contig_name][len(mod_rc[contig_name]) - pos] = ipd_ratio

    # Aggregate IPD ratio by motif sequence for each target object (bin or MGE) in parallel
    exe_pool = Pool(args.n_core)
    arg = [(t_name,
            list(contigs.keys()),
            list(contigs.values()),
            {c: mod_f[c] for c in contigs.keys()},
            {c: mod_rc[c] for c in contigs.keys()},
            motifs)
           for t_name, contigs in targets.items()]
    # NOTE: IPD vector = IPD ratio for each motif
    ipd_vectors = {t_name: ipd_vector
                   for t_name, ipd_vector in exe_pool.map(calc_ipdr, arg)}

    # Output as mean IPD ratio matrix of shape [target] x [motif]
    with open("ipdr_matrix.csv", 'w') as f:
        f.write('\t'.join([f"{seq} {pos} {met_type}" for seq, pos, met_type in motifs]))   # header
        f.write('\n')
        for t_name, ipd_vector in ipd_vectors.items():
            f.write('\t'.join([t_name.split('/')[-1],
                               '\t'.join(map(str, ipd_vector))]))
            f.write('\n')
