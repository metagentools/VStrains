#/bin/zsh
conda init zsh
conda activate spades-hapConstruction-env

# simulation
# 5 HIV
python eval_script/quast_evaluation.py -d benchmark/sim_result/5HIV/ -ref benchmark/strains/5-strain-HIV.fasta -o quast5hiv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py 

# 6 POLIO
python eval_script/quast_evaluation.py -d benchmark/sim_result/6POLIO/ -ref benchmark/strains/6-strain-polio.fasta -o quast6polio/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# 10 HCV
python eval_script/quast_evaluation.py -d benchmark/sim_result/10HCV/ -ref benchmark/strains/10-strain-HCV.fasta -o quast10hcv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# 15 ZIKV
python eval_script/quast_evaluation.py -d benchmark/sim_result/15ZIKV/ -ref benchmark/strains/15-strain-ZIKV.fasta -o quast15zikv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py 

# real
# 5 HIV labmix
python eval_script/quast_evaluation.py -d benchmark/real_result/5HIV_labmix/ -ref benchmark/5-virus-mix/data/REF.fasta -o quast5hiv_labmix/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# SARS-COV2 500x
python eval_script/quast_evaluation.py -d benchmark/real_result/SARS_COV2/500x/ -ref benchmark/real_result/SARS_COV2/Ref/ref_500x.fasta -o quast2sarscov_500x/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# SARS-COV2 4000x
python eval_script/quast_evaluation.py -d benchmark/real_result/SARS_COV2/4000x/ -ref benchmark/real_result/SARS_COV2/Ref/ref_4000x.fasta -o quast2sarscov_4000x/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

