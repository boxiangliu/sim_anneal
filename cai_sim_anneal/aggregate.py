import dill
from cai_sim_anneal.optimize import Plotter
import pandas as pd
import os
from utils.utils import cfg_file, load_config
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, nargs="+", help="input pickle files")
parser.add_argument("--fig", type=str, help="output figure file")
parser.add_argument("--organism", type=str,
                    help="human or yeast", default="human")
parser.add_argument("--protein", type=str,
                    help="spike or egfp", default="spike")
parser.add_argument("--cdsfold", type=str,
                    help="path to cdsfold result", default=None)
args = parser.parse_args()

# setup:
cfg = load_config(cfg_file)
if args.organism == "human":
    ref_p_file = cfg.DATA.RAW.REF_P.SPIKE.HUMAN

elif args.organism == "yeast":

    if args.protein == "spike":
        ref_p_file = cfg.DATA.RAW.REF_P.SPIKE.YEAST
    elif args.protein == "egfp":
        ref_p_file = cfg.DATA.RAW.REF_P.EGFP.YEAST
    else:
        raise ValueError(f"{args.protein} not implemented yet")

else:
    raise ValueError(f"{args.organism} not implemented yet")

ref_points = pd.read_csv(ref_p_file)
plotter = Plotter(ref_points)
color = ["orange", "green", "red", "purple", "brown", "pink",
         "gray", "olive", "cyan", "magenta", "black", "yellow", "teal",
         "chocolate", "salmon", "lightblue", "orchid", "indigo", 
         "tan", "lime", "navy", "gold", "lavender", "khaki"]


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

    lambda_ = extract_lambda_from_filename(pkl_file)
    plotter.add_points(sim_anneal, name=f"lambda: {lambda_}",
                       color=color[i], shape="^")

if args.cdsfold:
    cdsfold = pd.read_csv(args.cdsfold)
    plotter.add_points(cdsfold, name="CDSfold", color="blue", shape="s")

# make plot:
p = plotter.plot_cai_vs_mfe()
fig_file = os.path.join(
    cfg.DATA.PROCESSED.CAI_ANNEAL, args.fig)
p.save(fig_file)
