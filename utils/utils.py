import yaml
from easydict import EasyDict as edict
import subprocess
from collections import defaultdict
import numpy as np

cfg_file = "config.yaml"
def load_config(cfg_file):
    with open(cfg_file) as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
    return edict(cfg)

def read_fasta(file):
    seqs = {}
    with open(file, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip().replace(">","")
                id_ = header.split("|")[0]
            else:
                seq = line.strip()
                seqs[id_] = seq
    return seqs

def get_mfe(seq, folding_software):
    cmd = f"echo {seq} | {folding_software}"
    output = subprocess.run(cmd, capture_output=True, shell=True)
    mfe = output.stdout.decode("utf-8").split("\n")[1].split(" (")[1]
    return float(mfe.replace(")", ""))


class CAI(object):
    def __init__(self, codon_path):
        self.codon_path = codon_path
        self.nucs = set(["A", "U", "C", "G"])
        self.codon_table, self.max_aa_table = \
            self.parse_codon_table()

    def parse_codon_table(self):
        codon_table = defaultdict()
        max_aa_table = defaultdict(lambda: 0)
        for index, line in enumerate(open(self.codon_path).readlines()):
            if index == 0:
                continue

            codon, aa, fraction = line.strip().split(',')
            for e in codon:
                assert(e in self.nucs)
            codon_table[codon] = (float(fraction), aa)
            max_aa_table[aa] = max(max_aa_table[aa], float(fraction))

        return codon_table, max_aa_table

    def get_cai(self, rna_seq):
        assert(len(rna_seq) % 3 == 0)
        protein_length = len(rna_seq) // 3

        rna_seq = rna_seq.upper()
        rna_seq = rna_seq.replace("T","U")

        log_cai = 0.
        for index in range(protein_length):
            start = index * 3
            end = start + 3

            three = rna_seq[start:end]
            f_ci, aa = self.codon_table[three]

            f_c_max = self.max_aa_table[aa]
            w_i = f_ci / f_c_max
            log_cai += np.log(w_i)

        log_cai = log_cai / protein_length

        return np.exp(log_cai)
