# Random CAI optimization:
python cai_random/optimize.py "results_step-0.01_reps-50_org-human" --n_reps 50 --step_size 0.01 --organism human
python cai_random/optimize.py "results_step-0.01_reps-25_org-yeast" --n_reps 25 --step_size 0.01 --organism yeast


# CAI simulated annealing: 
python cai_sim_anneal/optimize.py 100 30 "results_it-100_factor-0.001_lamb-30" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0

python cai_sim_anneal/optimize.py 100 30 "results_it-100_factor-0.1_lamb-30" \
--objective "min" --factor 0.1 --anneal_schedule "linear" --seed 0

python cai_sim_anneal/optimize.py 100 30 "results_it-100_factor-0.001_lamb-30_linfold" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


python cai_sim_anneal/optimize.py 100 3 "results_it-100_factor-0.001_lamb-3_linfold" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


python cai_sim_anneal/optimize.py 1000 30 "results_it-1000_factor-0.001_lamb-30_linfold" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


python cai_sim_anneal/optimize.py 5000 30 "results_it-5000_factor-0.001_lamb-30_linfold" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


python cai_sim_anneal/optimize.py 1000 30 "results_it-1000_factor-0.001_lamb-30_linfold_onlyHigherCAI" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


python cai_sim_anneal/optimize.py 2000 15 "results_it-2000_factor-0.001_lamb-15_linfold_onlyHigherCAI" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


python cai_sim_anneal/optimize.py 2000 7 "results_it-2000_factor-0.001_lamb-15_linfold_onlyHigherCAI" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


python cai_sim_anneal/optimize.py 2000 3 "results_it-2000_factor-0.001_lamb-15_linfold_onlyHigherCAI" \
--objective "min" --factor 0.001 --anneal_schedule "linear" --seed 0 --folding "LinearFold"


for lambda in 0.5 1 2 3 4 6 8 10; do
    sbatch -p TitanXx8_slong,M40x8_slong,1080Ti_slong -o lambda-${lambda}.log --wrap "python cai_sim_anneal/optimize.py 30000 $lambda results_it-30000_factor-0.001_lamb-${lambda}_linfold_onlyHigherCAI --objective min --factor 0.001 --anneal_schedule linear --seed 0 --folding LinearFold"
done

for lambda in 3.5 4.5 5 5.5 7; do
    sbatch -p TitanXx8_slong,M40x8_slong,1080Ti_slong -o lambda-${lambda}.log --wrap "python cai_sim_anneal/optimize.py 30000 $lambda results_it-30000_factor-0.001_lamb-${lambda}_linfold_onlyHigherCAI --objective min --factor 0.001 --anneal_schedule linear --seed 0 --folding LinearFold"
done


for lambda in 0.5 1 2 3 3.5 4 4.5 5 5.5 6 7 8 10; do
    sbatch -p TitanXx8_slong,M40x8_slong,1080Ti_slong -o lambda-${lambda}.log --wrap "python cai_sim_anneal/optimize.py 30000 $lambda results_it-30000_factor-0.001_lamb-${lambda}_linfold_onlyHigherCAI_yeast --objective min --factor 0.001 --anneal_schedule linear --seed 0 --folding LinearFold --organism yeast"
done

python cai_sim_anneal/aggregate.py --input ../data/processed/cai_sim_anneal/results_it-30000_factor-0.001_lamb-{0.5,1,2,3,3.5,4,4.5,5,5.5,6,7,8,10}_linfold_onlyHigherCAI.log --fig test.pdf
