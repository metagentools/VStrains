#/bin/zsh
conda init zsh
conda activate spades-hapConstruction-env

python eval_script/quast_evaluation.py -d benchmark/sim_result/5HIV/ -ref benchmark/strains/5-strain-HIV.fasta -o quast5hiv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py 
# # 5 HIV
python src/__init__.py -a spades -g benchmark/fastq/5-strain-HIV-20000x/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/5-strain-HIV-20000x/output_careful/contigs.paths -o acc_5_hiv_careful/ -d
# # python eval_script/quast_evaluation.py -cs benchmark/fastq/5-strain-HIV-20000x/output_careful/contigs.fasta ~/Desktop/benchmark/shortread/vgflow+savage/vgflow5hiv/haps.final.fasta acc_5_hiv_careful/strain.fasta -ref benchmark/strains/5-strain-HIV.fasta -o quast5hiv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py 

# # 6 POLIO
# python src/__init__.py -a spades -g benchmark/fastq/6-strain-poliovirus/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/6-strain-poliovirus/output_careful/contigs.paths -o acc_6_polio/ -d
# # python eval_script/quast_evaluation.py -cs benchmark/fastq/6-strain-poliovirus/output_careful/contigs.fasta ~/Desktop/benchmark/shortread/vgflow+savage/vgflow6polio/haps.final.fasta acc_6_polio/strain.fasta -ref benchmark/strains/6-strain-polio.fasta -o quast6polio/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py 

# # 10 HCV
# python src/__init__.py -a spades -g benchmark/fastq/10-strain-HCV-20000x/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/10-strain-HCV-20000x/output_careful/contigs.paths -o acc_10_hcv/ -d
# # python eval_script/quast_evaluation.py -cs benchmark/fastq/10-strain-HCV-20000x/output_careful/contigs.fasta ~/Desktop/benchmark/shortread/vgflow+savage/vgflow10hcv/haps.final.fasta acc_10_hcv/strain.fasta -ref benchmark/strains/10-strain-HCV.fasta -o quast10hcv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py 

# # 15 ZIKV
# python src/__init__.py -a spades -g benchmark/fastq/15-strain-ZIKV-20000x/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/15-strain-ZIKV-20000x/output_careful/contigs.paths -o acc_15_zikv_careful/ -d
# # python eval_script/quast_evaluation.py -cs benchmark/fastq/15-strain-ZIKV-20000x/output_careful/contigs.fasta ~/Desktop/benchmark/shortread/vgflow+savage/vgflow15zikv/haps.final.fasta acc_15_zikv_careful/strain.fasta -ref benchmark/strains/15-strain-ZIKV.fasta -o quast15zikv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py 

# python src/__init__.py -a spades -g ~/Desktop/5-virus-mix/output_5hiv_labmix/assembly_graph_after_simplification.gfa -p ~/Desktop/5-virus-mix/output_5hiv_labmix/contigs.paths -o acc_5hiv_labmix/ -d
# # python eval_script/quast_evaluation.py -cs ~/Desktop/5-virus-mix/output_5hiv_labmix/contigs.fasta ../vg-flow/vgflow5hiv_labmix_savage/contigs_stage_c.fasta acc_5hiv_labmix/strain.fasta ../vg-flow/vgflow5hiv_labmix_spades/haps.final.fasta ../vg-flow/vgflow5hiv_labmix_savage/haps.final.fasta -ref ~/Desktop/5-virus-mix/data/REF.fasta -o quast_5hiv_labmix/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

