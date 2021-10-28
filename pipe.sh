# Random CAI optimization:
python cai_random/mutate.py


# CAI simulated annealing: 
python cai_sim_anneal/optimize.py 100 30 "results_it-100_factor-0.001_lamb-30" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0

python cai_sim_anneal/optimize.py 100 30 "results_it-100_factor-0.1_lamb-30" \
--objective "min" --factor 0.1 --anneal_schedule "linear" --seed 0

python cai_sim_anneal/optimize.py 100 30 "results_it-100_factor-0.1_lamb-30_linfold" \
--objective "min" --factor 0.1 --anneal_schedule "linear" --seed 0 --folding "LinearFold"