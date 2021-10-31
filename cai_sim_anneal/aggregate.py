import pickle as pkl
from cai_sim_anneal.optimize import Plotter

ref_points = pd.read_csv(ref_p_file)
plotter = Plotter(ref_points)

color = ["orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan", "magenta", "black", "yellow"]

for i, lambda_ in enumerate([0.5, 1, 2, 3, 4, 6, 8, 10]):
    pkl_file = os.path.join(cfg.DATA.PROCESSED.CAI_ANNEAL, f"results_it-30000_factor-0.001_lamb-{lambda_}_linfold_onlyHigherCAI.pkl")
    with open(pkl_file, "rb") as f:
        results = pkl.load(f)["results"]

    sim_anneal = [(k, v["MFE"], v["CAI"]) for k, v in results.items()]
    sim_anneal = pd.DataFrame(
        sim_anneal, columns=["Iteration", "MFE", "CAI"])

    plotter.add_points(sim_anneal, name="Simulated Annealing",
                       color=color[i], shape="^")


# make plot
p = plotter.plot_cai_vs_mfe()
fig_file = os.path.join(
    cfg.DATA.PROCESSED.CAI_ANNEAL, "test" + ".pdf")