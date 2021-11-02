import dill
from cai_sim_anneal.optimize import Plotter
import pandas as pd
import os
from utils.utils import cfg_file, load_config
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, nargs="+", help="input pickle files")
parser.add_argument("--fig", type=str, help="output figure file")
args = parser.parse_args()

# setup:
cfg = load_config(cfg_file)
ref_p_file = cfg.DATA.RAW.REF_P.HUMAN
ref_points = pd.read_csv(ref_p_file)
plotter = Plotter(ref_points)
color = ["orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan", "magenta", "black", "yellow"]


def extract_lambda_from_filename(filename):
    split_filename = filename.strip().split("_")
    for i in split_filename:
        if i.startswith("lamb"):
            return i.replace("lamb-", "")

# load different lambdas:
for i, pkl_file in enumerate(args.input):
    with open(pkl_file, "rb") as f:
        results = dill.load(f)["results"]

    sim_anneal = [(k, v["MFE"], v["CAI"]) for k, v in results.items()]
    sim_anneal = pd.DataFrame(
        sim_anneal, columns=["Iteration", "MFE", "CAI"])

    lambda_ = extract_lambda_from_file(pkl_file)
    plotter.add_points(sim_anneal, name=f"lambda: {lambda_}",
                       color=color[i], shape="^")


# make plot:
p = plotter.plot_cai_vs_mfe()
fig_file = os.path.join(
    cfg.DATA.PROCESSED.CAI_ANNEAL, args.fig)
p.save(fig_file)