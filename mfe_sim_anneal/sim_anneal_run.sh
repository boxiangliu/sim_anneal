sim_anneal(){
	python3 sim_anneal.py \
	--protein membrane.protein \
	--length $1 \
	--out_dir $2 \
	-c $3 --seed $4
}

export -f sim_anneal

parallel -j 4 sim_anneal 223 runs/length-223/ {1} {2} ::: 1 10 100 1000 ::: {1..3}
parallel -j 4 sim_anneal 70 runs/length-70/ 1 {1} ::: {1..100}
parallel -j 4 sim_anneal 100 runs/length-100/ 1 {1} ::: {1..100}