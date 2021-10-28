# Random CAI optimization:
python cai_random/mutate.py


# CAI simulated annealing: 
python cai_sim_anneal/optimize.py 100 0.5 "results_it-100_lamb-0.5" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0