# Random CAI optimization:
python cai_random/optimize.py "results_step-0.01_reps-50_org-human" --n_reps 50 --step_size 0.01 --organism human
python cai_random/optimize.py "results_step-0.01_reps-25_org-yeast" --n_reps 25 --step_size 0.01 --organism yeast


###########################
# CAI simulated annealing:
###########################
# Human
for lambda in 0.5 1 2 3 3.5 4 4.5 5 5.5 6 7 8 10; do
    sbatch -p TitanXx8_slong,M40x8_slong,1080Ti_slong -o lambda-${lambda}.log --wrap "python cai_sim_anneal/optimize.py 30000 $lambda results_it-30000_factor-0.001_lamb-${lambda}_linfold_onlyHigherCAI --objective min --factor 0.001 --anneal_schedule linear --seed 0 --folding LinearFold"
done

python cai_sim_anneal/aggregate.py --input ../data/processed/cai_sim_anneal/results_it-30000_factor-0.001_lamb-{0.5,1,2,3,3.5,4,4.5,5,5.5,6,7,8,10}_linfold_onlyHigherCAI.pkl \
    --fig aggregate_human.pdf --organism human

# Yeast
for lambda in 0.5 1 2 3 3.5 4 4.5 5 5.5 6 7 8 10; do
    sbatch -p TitanXx8_slong,M40x8_slong,1080Ti_slong -o lambda-${lambda}.log --wrap "python cai_sim_anneal/optimize.py 30000 $lambda results_it-30000_factor-0.001_lamb-${lambda}_linfold_onlyHigherCAI_yeast --objective min --factor 0.001 --anneal_schedule linear --seed 0 --folding LinearFold --organism yeast"
done

python cai_sim_anneal/aggregate.py --input ../data/processed/cai_sim_anneal/results_it-30000_factor-0.001_lamb-{0.5,1,2,3,3.5,4,4.5,5,5.5,6,7,8,10}_linfold_onlyHigherCAI_yeast.pkl \
    --fig aggregate_yeast.pdf --organism yeast
