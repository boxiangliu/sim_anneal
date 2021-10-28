from collections import defaultdict
import subprocess
import os
import numpy as np
from tqdm import tqdm
from utils.utils import read_fasta, cfg_file, load_config, CAI, \
    read_coding_wheel, get_equivalent_codons
import pickle as pkl
import dill
import pandas as pd
from plotnine import ggplot, geom_point, geom_line, aes, theme_bw, \
    scale_color_manual, scale_shape_manual, theme, element_blank
import argparse
import logging
import time

parser = argparse.ArgumentParser()
parser.add_argument("iteration", type=int, help="number of iterations to run")
parser.add_argument("lambd", type=float, help="CAI coefficient")
parser.add_argument("out_file", type=str, help="path to output file")
parser.add_argument("--objective", type=str,
                    help="either min or max", default="min")
parser.add_argument("--factor", type=float,
                    help="scaling factor for temperature", default=0.001)
parser.add_argument("--anneal_schedule", type=str,
                    help="either linear or geometric", default="linear")
parser.add_argument("--alpha", type=float,
                    help="scale factor for temperature update", default=None)
parser.add_argument("--seed", type=int, help="random seed", default=None)
args = parser.parse_args()


class RNA():

    def __init__(self, seq, equi_codons, codon_freq, folding_cmd):
        self.equi_codons = equi_codons
        self.codon_freq = codon_freq
        self.folding_cmd = folding_cmd

        self.set_rna(seq)
        self.cai_calc = CAI(codon_freq)

        self.nt_len = len(seq)
        assert self.nt_len % 3 == 0
        self.pt_len = len(seq) // 3

    def set_rna(self, seq):
        seq = seq.upper()
        seq = seq.replace("T", "U")
        self.rna = seq

    def get_mfe(self):
        cmd = f"echo {self.rna} | {self.folding_cmd}"
        output = subprocess.run(cmd, capture_output=True, shell=True)
        mfe = output.stdout.decode("utf-8").split("\n")[1].split(" (")[1]
        return float(mfe.replace(")", ""))

    def get_cai(self):
        return self.cai_calc.get_cai(self.rna)

    def get_score(self, lambda_):
        mfe = self.get_mfe()
        cai = self.get_cai()
        score = mfe - self.pt_len * lambda_ * np.log(cai)
        return mfe, cai, score

    def mutate(self):

        # get equivalent codons:
        codon_set = {}
        while len(codon_set) == 0:
            mut_id = np.random.choice(self.pt_len)  # choose a random codon
            start = mut_id * 3
            end = start + 3
            codon = self.rna[start:end]
            codon_set = self.equi_codons[codon]

        # mutate:
        mutation = np.random.choice(list(codon_set))
        mut_seq = list(self.rna)
        mut_seq[start:end] = mutation
        # print(f"{codon} -> {mutation}")
        self.rna = "".join(mut_seq)


class SimAnnealer(object):

    def __init__(self, model, iteration, lambda_, objective="max", factor=1,
                 anneal_schedule="linear", alpha=None, seed=None):
        self.model = model
        self.iteration = iteration
        self.lambda_ = lambda_
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
            if self.anneal_schedule == "linear":
                self.alpha = 1 / self.iteration
            else:
                raise ValueError(
                    "alpha must be set when objective is not linear")

    def update_temperature(self, T):
        if self.anneal_schedule == "linear":
            T -= self.alpha
        elif self.anneal_schedule == "geometric":
            T *= self.alpha
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

        start = time.time()

        mfe, cai, old_score = self.model.get_score(self.lambda_)
        old_seq = self.model.rna
        T = 1
        for i in tqdm(range(self.iteration)):

            self.model.mutate()
            mfe, cai, new_score = self.model.get_score(self.lambda_)

            end = time.time()
            elapsed = end - start
            start = end

            if self.better(old_score, new_score) or \
                    (np.random.uniform() <= self.get_prob(old_score, new_score, T)):
                old_score = new_score
                old_seq = self.model.rna
                self.results[i] = {"seq": self.model.rna,
                                   "score": new_score, "CAI": cai,
                                   "MFE": mfe, "lambda": self.lambda_}
                logging.info(f"ELAPSED: {elapsed}\tSCORE: {new_score}\tACCEPT: YES\tCAI: {cai}\tMFE: {mfe}")
            else:
                self.model.rna = old_seq
                logging.info(f"ELAPSED: {elapsed}\tSCORE: {new_score}\tACCEPT: NO\tCAI: {cai}\tMFE: {mfe}")

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

        with open(file + ".pkl", "wb") as f:
            dill.dump(results, f)


class Plotter(object):

    def __init__(self, ref_points):
        self.points = {}
        self.points["Reference"] = (ref_points, "blue", "o")

    def add_points(self, points, name, color, shape):
        self.points[name] = (points, color, shape)

    def plot_cai_vs_iteration(self):
        p = (ggplot(data=self.points["Simulated Annealing"][0],
                    mapping=aes(x="Iteration", y="CAI"))
             + geom_point()
             + theme_bw())
        return p

    def plot_cai_vs_mfe(self):
        plot_df = defaultdict(list)
        color_map = {}
        shape_map = {}
        for name in self.points:
            (points, color, shape) = self.points[name]

            plot_df["MFE"] += points["MFE"].tolist()
            plot_df["CAI"] += points["CAI"].tolist()
            plot_df["name"] += [name] * points.shape[0]
            color_map[name] = color
            shape_map[name] = shape
        plot_df = pd.DataFrame(plot_df)

        p = (ggplot(data=plot_df, mapping=aes(x="MFE", y="CAI", color="name", shape="name"))
             + geom_point()
             + theme_bw()
             + scale_color_manual(values=color_map)
             + scale_shape_manual(values=shape_map)
             + theme(legend_position="top", legend_title=element_blank()))

        return p


def run(args, cfg):
    iteration = args.iteration
    lambda_ = args.lambd
    out_file = args.out_file
    objective = args.objective
    factor = args.factor
    anneal_schedule = args.anneal_schedule
    alpha = args.alpha
    seed = args.seed

    seqs = read_fasta(cfg.DATA.RAW.SPIKE)
    mfe_seq = seqs["lambda_0"]

    codon_table = read_coding_wheel(cfg.DATA.RAW.CODON_TABLE)
    equi_codons = get_equivalent_codons(codon_table)
    model = RNA(mfe_seq, equi_codons, cfg.DATA.RAW.CODON_FREQ,
                folding_cmd=cfg.BIN.RNAFOLD)

    annealer = SimAnnealer(model, iteration=iteration, lambda_=lambda_, objective=objective,
                           factor=factor, anneal_schedule=anneal_schedule, seed=seed)
    print(annealer)

    annealer.anneal()
    annealer.save(os.path.join(cfg.DATA.PROCESSED.CAI_ANNEAL, out_file))


def main():
    # iteration = 50000
    # lambda_ = 0.5
    # objective = "min"
    # factor = 0.001
    # anneal_schedule = "linear"
    # alpha = None
    # seed = 0
    # out_file = "results.pkl"
    cfg = load_config(cfg_file)
    log_file = os.path.join(cfg.DATA.PROCESSED.CAI_ANNEAL, args.out_file + ".log")
    logging.basicConfig(
        level=logging.INFO,
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()])
    run(args, cfg)

if __name__ == "__main__":
    main()
