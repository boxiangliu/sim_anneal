from collections import defaultdict
import subprocess
import os
import numpy as np
from tqdm import tqdm
from utils.utils import read_fasta, cfg_file, load_config, CAI


class RNA():

    def __init__(self, seq, codon_table, codon_freq):
        self.codon_table = codon_table
        self.codon_freq = codon_freq
        self.rna = seq
        self.cai_calc = CAI(codon_freq)

    def get_mfe(self, folding_software):
        cmd = f"echo {self.rna} | {folding_software}"
        output = subprocess.run(cmd, capture_output=True, shell=True)
        mfe = output.stdout.decode("utf-8").split("\n")[1].split(" (")[1]
        return float(mfe.replace(")", ""))

    def get_cai(self):
        return self.cai_calc.get_cai(self.rna)

    def get_score(self):
        return self.get_cai()

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

    def __init__(self, model, iteration, objective="max", factor=1000,
                 anneal_schedule="linear", alpha=None, seed=None):
        self.model = model
        self.iteration = iteration
        self.objective = objective
        self.factor = factor
        self.anneal_schedule = anneal_schedule
        self.set_alpha(alpha)
        self.seed = seed

        self.results = dict()

    def __repr__(self):
        attr = [f"Model:             {self.model}",
                f"Iteration:         {self.iteration}",
                f"Objective:         {self.objective}",
                f"Factor:            {self.factor}",
                f"Anneal Schedule:   {self.anneal_schedule}",
                f"Alpha:             {self.alpha}",
                f"Seed:              {self.seed}"]
        return "\n".join(attr)

    def set_alpha(self, alpha):
        if alpha != None:
            self.alpha = alpha
        else:
            if anneal_schedule == "linear":
                self.alpha = 1 / self.iteration
            else:
                raise ValueError(
                    "alpha must be set when objective is not linear")

    def update_temperature(self, T):
        if self.anneal_schedule == "linear":
            T -= alpha
        elif self.anneal_schedule == "geometric":
            T *= alpha
        else:
            raise ValueError(f"{method} not yet implemented!")
        return T

    def better(self, old_score, new_score):
        objective = self.objective

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

    def get_prob(self, old_score, new_score, T):
        return np.exp(-abs(old_score - new_score) / (self.factor * T))

    def anneal(self):
        if self.seed:
            np.random.seed(self.seed)

        old_score = self.model.get_score()
        old_seq = model.rna
        T = 1
        for i in tqdm(range(self.iteration)):

            model.mutate()
            new_score = model.get_score()

            if self.better(old_score, new_score) or \
                    (np.random.uniform() <= self.get_prob(old_score, new_score, T)):

                old_score = new_score
                old_seq = model.rna
                self.results[(proportion, seed)] = [mut_seq, mfe, cai]

            else:
                model.rna = old_seq

            T = self.update_temperature(T)

    def save(self, file):
        results = {"model": self.model,
                   "iteration": self.iteration,
                   "objective": self.objective,
                   "factor": self.factor,
                   "anneal_schedule": self.anneal_schedule,
                   "alpha": self.alpha,
                   "seed": self.seed,
                   "results": self.results}

        with open(file, "wb") as f:
            pkl.dump(results, f)


def run(iteration, objective, factor, anneal_schedule, alpha, seed):
    cfg = load_config(cfg_file)
    seqs = read_fasta(cfg.DATA.RAW.SPIKE)
    mfe_seq = seqs["lambda_0"]
    model = RNA(mfe_seq, cfg.DATA.RAW.CODON_TABLE,
                cfg.DATA.RAW.CODON_FREQ)

    annealer = SimAnnealer(model, iteration=1000, factor=1000, seed=0)
    print(annealer)

    annealer.anneal()
    annealer.save(cfg.DATA.PROCESSED.CAI_ANNEAL)


def main():
    iteration = 1000
    objective = "max"
    factor = 1000
    anneal_schedule = "linear"
    alpha = None
    seed = 0
    run(iteration, objective, factor, anneal_schedule, alpha, seed)

if __name__ == "__main__":
    main()
