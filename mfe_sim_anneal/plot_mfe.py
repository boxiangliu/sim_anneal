import glob
import os
from collections import defaultdict
import pandas as pd
import seaborn as sns
import click

def get_dataframe(in_dir):
	proto_df = defaultdict(list)
	for subdir in glob.glob(f"{in_dir}/*"):
		print(subdir)
		c, seed = os.path.basename(subdir).split("_")
		c = int(c.split("-")[1])
		seed = int(seed.split("-")[1])

		mfe_fn = f"{subdir}/mfe.txt"

		with open(mfe_fn, "r") as f:
			for line in f:
				iteration, mfe = line.strip().split("\t")
				proto_df["iteration"].append(int(iteration))
				proto_df["mfe"].append(float(mfe))
				proto_df["c"].append(c)
				proto_df["seed"].append(seed)
	df = pd.DataFrame(proto_df)
	return df

def select_mfe(df):
	min_mfe = df["mfe"].min()
	seed = df[df["mfe"] == min_mfe]["seed"].unique()[0]
	return df[df["seed"] == seed]

def p1(df, out_dir):
	c = df["c"].unique()[0]
	seed = df["seed"].unique()[0]
	mfe = df["mfe"].min()
	fig = sns.lineplot(x="iteration", y="mfe", \
		data=df).set_title(f"c = {c}; seed = {seed}; mfe = {mfe}")
	fig.get_figure().savefig(f"{out_dir}/c-{c}_seed-{seed}.png")
	fig.get_figure().clf()


@click.command()
@click.option("-i", "--in_dir", type=str, help="Input directory", \
	required=True)
@click.option("-o", "--out_dir", type=str, help="Output directory", \
	required=True)
def main(in_dir, out_dir):
	os.makedirs(out_dir, exist_ok=True)

	df = get_dataframe(in_dir)
	df = select_mfe(df)
	p1(df, out_dir)


if __name__ == "__main__":
	main()