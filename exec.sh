#/bin/zsh
conda init zsh
conda activate spades-hapConstruction-env

# # 5 HIV
python src/vsAware/__init__.py -a spades -g benchmark/fastq/5-strain-HIV-20000x/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/5-strain-HIV-20000x/output_careful/contigs.paths -o acc_5_hiv/ -fwd benchmark/fastq/5-strain-HIV-20000x/forward.fastq -rve benchmark/fastq/5-strain-HIV-20000x/reverse.fastq
# python src/vsAware/utils/quast_evaluation.py -d benchmark/sim_result/5HIV/ -ref benchmark/strains/5-strain-HIV.fasta -o quast5hiv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# # 6 POLIO
# python src/vsAware/__init__.py -a spades -g benchmark/fastq/6-strain-poliovirus/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/6-strain-poliovirus/output_careful/contigs.paths -o acc_6_polio/ -e1 6polio_pe_info_relax/pe_info -e2 6polio_pe_info_relax/st_info -fwd benchmark/fastq/6-strain-poliovirus/all.forward.fastq -rve benchmark/fastq/6-strain-poliovirus/all.reverse.fastq
# python src/vsAware/utils/quast_evaluation.py -d benchmark/sim_result/6POLIO/ -ref benchmark/strains/6-strain-polio.fasta -o quast6polio/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# # 10 HCV
# python src/vsAware/__init__.py -a spades -g benchmark/fastq/10-strain-HCV-20000x/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/10-strain-HCV-20000x/output_careful/contigs.paths -o acc_10_hcv/ -e1 10hcv_pe_info_relax/pe_info -e2 10hcv_pe_info_relax/st_info -fwd benchmark/fastq/10-strain-HCV-20000x/forward.fastq -rve benchmark/fastq/10-strain-HCV-20000x/reverse.fastq
# python src/vsAware/utils/quast_evaluation.py -d benchmark/sim_result/10HCV/ -ref benchmark/strains/10-strain-HCV.fasta -o quast10hcv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# # 15 ZIKV
# python src/vsAware/__init__.py -a spades -g benchmark/fastq/15-strain-ZIKV-20000x/output_careful/assembly_graph_after_simplification.gfa -p benchmark/fastq/15-strain-ZIKV-20000x/output_careful/contigs.paths -o acc_15_zikv/ -e1 15zikv_pe_info_relax/pe_info -e2 15zikv_pe_info_relax/st_info -fwd benchmark/fastq/15-strain-ZIKV-20000x/forward.fastq -rve benchmark/fastq/15-strain-ZIKV-20000x/reverse.fastq  
# python src/vsAware/utils/quast_evaluation.py -d benchmark/sim_result/15ZIKV/ -ref benchmark/strains/15-strain-ZIKV.fasta -o quast15zikv/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py

# # 5 HIV-Labmix
# python src/vsAware/__init__.py -a spades -g benchmark/5-virus-mix/output_5hiv_labmix/assembly_graph_after_simplification.gfa -p benchmark/5-virus-mix/output_5hiv_labmix/contigs.paths -o acc_5hiv_labmix/ -e1 5hiv_labmix_pe_info_relax/pe_info -e2 5hiv_labmix_pe_info_relax/st_info -fwd benchmark/fastq/real-5-strain-HIV-20000x/SRR961514_1.fp.fastq -rve benchmark/fastq/real-5-strain-HIV-20000x/SRR961514_2.fp.fastq
# python src/vsAware/utils/quast_evaluation.py -d benchmark/real_result/5HIV_labmix/ -ref benchmark/5-virus-mix/data/REF.fasta -o quast5hiv_labmix/ -quast ~/bio_tools/quast-5.1.0rc1/metaquast.py
