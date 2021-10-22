from collections import defaultdict
import subprocess
import os
import math
import random
import ipdb
import click
import numpy as np

CODING_WHEEL_FN = "coding_wheel.txt"
PROTEIN_FN = "membrane.protein"


def read_coding_wheel(fn):
    codon_table = defaultdict(list)
    with open(fn, "r") as f:
        for line in f:
            split_line = line.strip().split("\t")
            aa = split_line[0]
            for codon in split_line[1:]:
                first, second, thirds = codon.split(" ")
                for third in thirds:
                    codon_table[aa].append(first + second + third)
    return codon_table


def read_protein(fn):
    protein = []
    with open(fn, "r") as f:
        line = f.readline()
        for aa in line.strip().split(" "):
            protein.append(aa)
    return protein


class RNA():

    def __init__(self, seq, codon_table):
        self.codon_table = codon_table
        self.rna = seq

    def get_mfe(self, folding_software):
        cmd = f"echo {self.rna} | {folding_software}"
        output = subprocess.run(cmd, capture_output=True, shell=True)
        mfe = output.stdout.decode("utf-8").split("\n")[1].split(" (")[1]
        return float(mfe.replace(")", ""))

    def get_cai(self):
        pass

    def mutate(self):
        indices = [x for x in range(len(self.protein))]
        codon_set = [""]
        while len(codon_set) == 1:
            index = random.choice(indices)
            aa = self.protein[index]
            codon_set = self.codon_table[aa]

        codon = self.rna[index]
        mutables = [x for x in codon_set if x != codon]
        mutation = random.choice(mutables)
        self.rna[index] = mutation
        print(f"{codon} -> {mutation}")


class SimAnnealer(object):

    def __init__(self, model, objective):
        self.model = model
        self.results = dict()
        self.objective = objective

    def update_temperature(self, T, alpha, method="linear"):
        if method == "linear":
            T -= alpha
        elif method == "geometric":
            T *= alpha
        else:
            raise ValueError(f"{method} not yet implemented!")
        return T

    def better(self, old_score, new_score, objective="max"):
        if objective == "max":
            if new_score >= old_score:
                return True
            else:
                return False
        elif objective == "min":
            if new_score <= old_score:
                return True
            else:
                return False
        else:
            raise ValueError(f"{objective} not yet implemented!")

    def get_prob(self, old_score, new_score, c, T):
        return np.exp(-abs(old_score - new_score) / (c * T))

    def anneal(self, iteration, c, alpha, update_method="linear", seed=None):
        if seed:
            np.random.seed(seed)

        old_score = model.get_score()
        old_seq = model.rna
        T = 1
        for i in range(iteration):

            T = self.update_temperature(T, alpha, update_method)
            model.mutate()
            new_score = model.get_score()

            if self.better(score, new_score, objective=self.objective) or \
                    (np.random.uniform() <= self.get_prob(old_score, new_score, c, T)):

                old_score = new_score
                old_seq = model.rna
                self.results[(proportion, seed)] = [mut_seq, mfe, cai]

            else:
                model.rna = old_seq


def sim_anneal(protein, out_dir, iteration=1000, c=1000, seed=42):
    random.seed(seed)

    codon_table = read_coding_wheel(CODING_WHEEL_FN)
    rna = RNA(protein, codon_table)

    out_dir = f"{out_dir}/c-{c}_seed-{seed}/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    min_mfe = 0
    mfe_0 = 0
    rna_0 = ""
    f1 = open(f"{out_dir}/mfe.txt", "w")

    for i in range(iteration):
        T = 1 - i / iteration
        rna.mutate()
        mfe_1 = rna.get_mfe(f"{out_dir}/iter_{i:04}.fa")

        if mfe_0 >= mfe_1:
            print(f"[Iteration {i}]: MFE down; {mfe_0} -> {mfe_1}.")
            mfe_0 = mfe_1
            rna_0 = rna.rna
        else:
            p_update = math.exp((mfe_0 - mfe_1) / (c * T))

            if random.uniform(0, 1) <= p_update:
                print(f"[Iteration {i}]: MFE up; {mfe_0} -> {mfe_1}.")
                mfe_0 = mfe_1
                rna_0 = rna.rna
            else:
                rna.rna = rna_0
                rna.mfe = mfe_0
                print(f"[Iteration {i}]: MFE same.")

        if rna.mfe < min_mfe:
            min_mfe = rna.mfe
            with open(f"{out_dir}/min_mfe.fa", "w") as f2:
                f2.write(f">c-{c}_seed-{seed}_iteration-{i}\n")
                f2.write("".join(rna.rna) + "\n")

        f1.write(f"{i:04}\t{rna.mfe}" + "\n")

    f1.close()


@click.command()
@click.option("-p", "--protein", type=str, help="Protein file",
              required=True)
@click.option("-l", "--length", type=int, help="Select certain protein length.")
@click.option("-o", "--out_dir", type=str, help="Output directory.")
@click.option("-c", type=int, help="Temperature scaling parameter.",
              default=1, show_default=True)
@click.option("-s", "--seed", type=int, help="Random seed.", default=42,
              show_default=True)
@click.option("-i", "--iteration", type=int, help="Number of iterations.",
              default=1000, show_default=True)
def main(protein, length, out_dir, iteration, c, seed,):
    protein = read_protein(PROTEIN_FN)[:length]
    sim_anneal(protein, out_dir, iteration, c, seed)


if __name__ == "__main__":
    main()
