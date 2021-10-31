from utils.utils import read_fasta, cfg_file, load_config, CAI
import numpy as np
from tqdm import tqdm
import pickle as pkl
import os
import subprocess
from collections import defaultdict
import pandas as pd
from plotnine import ggplot, geom_point, geom_line, aes, theme_bw, \
    scale_color_manual, scale_shape_manual, theme, element_blank
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("out_file", type=str, help="output file")
parser.add_argument("--n_reps", type=int,
                    help="number of replications", default=50)
parser.add_argument("--step_size", type=float,
                    help="step size for proportion", default=0.01)
parser.add_argument("--organism", type=str,
                    help="organism: yeast or human", default="human")
args = parser.parse_args()


class CAIOptimizer(object):

    def __init__(self, mfe_seq, cai_seq):
        assert len(mfe_seq) % 3 == 0
        assert len(mfe_seq) == len(cai_seq)
        self.mfe_seq = mfe_seq
        self.cai_seq = cai_seq
        self.nt_len = len(self.mfe_seq)
        self.pt_len = self.nt_len // 3

    def mutate(self, proportion, seed=None):
        assert (proportion >= 0) and (proportion <= 1)
        num_mut = int(np.floor(self.pt_len * proportion))

        if seed:
            np.random.seed(seed)
        if num_mut == 0:
            return self.mfe_seq

        mut_ids = np.random.choice(self.pt_len, num_mut, replace=False)
        mut_seq = list(self.mfe_seq)
        for mut_id in mut_ids:
            start = 3 * mut_id
            mut_seq[start:start + 3] = self.cai_seq[start:start + 3]
        return "".join(mut_seq)

    def get_mfe(self, seq, folding_cmd):
        cmd = f"echo {seq} | {folding_cmd}"
        output = subprocess.run(cmd, capture_output=True, shell=True)
        mfe = output.stdout.decode("utf-8").split("\n")[1].split(" (")[1]
        return float(mfe.replace(")", ""))


class Plotter(object):

    def __init__(self, ref_points):
        self.points = {}
        self.points["Reference"] = (ref_points, "blue", "o")

    def add_points(self, points, name, color, shape):
        self.points[name] = (points, color, shape)

    def plot(self):
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
    # initialize:
    step_size = args.step_size
    n_reps = args.n_reps
    out_file = args.out_file
    organism = args.organism

    if organism == "human":
        codon_freq_file = cfg.DATA.RAW.CODON_FREQ.HUMAN
        spike_protein_file = cfg.DATA.RAW.SPIKE.HUMAN
        ref_p_file = cfg.DATA.RAW.REF_P.HUMAN
    elif organism == "yeast":
        codon_freq_file = cfg.DATA.RAW.CODON_FREQ.YEAST
        spike_protein_file = cfg.DATA.RAW.SPIKE.YEAST
        ref_p_file = cfg.DATA.RAW.REF_P.YEAST
    else:
        raise ValueError(f"{organism} not implemented!")

    seqs = read_fasta(spike_protein_file)
    mfe_seq = seqs["lambda_0"]
    cai_seq = seqs["lambda_inf"]

    cai_calc = CAI(codon_freq_file)

    # mutate sequence and calculate MFE & CAI:
    optimizer = CAIOptimizer(mfe_seq, cai_seq)
    results = dict()
    for proportion in tqdm(np.arange(0, 1.001, step_size)):
        for seed in np.arange(n_reps):
            mut_seq = optimizer.mutate(proportion, seed=seed)
            mfe = optimizer.get_mfe(seq=mut_seq, folding_cmd=cfg.BIN.RNAFOLD)
            cai = cai_calc.get_cai(mut_seq)
            results[(proportion, seed)] = [mut_seq, mfe, cai]

    pkl_file = os.path.join(
        cfg.DATA.PROCESSED.CAI_RANDOM, out_file + ".pkl")
    with open(pkl_file, "wb") as f:
        pkl.dump(results, f)

    # make MFE-CAI plot:
    ref_points = pd.read_csv(ref_p_file)
    cai_random = [(k[0], k[1], v[1], v[2]) for k, v in results.items()]
    cai_random = pd.DataFrame(
        cai_random, columns=["proportion", "seed", "MFE", "CAI"])

    plotter = Plotter(ref_points)
    plotter.add_points(cai_random, name="Random Optimization",
                       color="red", shape="^")
    p = plotter.plot()
    p.save(os.path.join(cfg.DATA.PROCESSED.CAI_RANDOM, out_file + ".pdf"))


def main():
    cfg = load_config(cfg_file)
    run(args, cfg)


if __name__ == "__main__":
    main()
