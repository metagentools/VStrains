#/bin/zsh
conda init zsh
conda activate spades-hapConstruction-env

# 5 HIV
time python src/hap_construction.py -gfa benchmark/fastq/5-strain-HIV-20000x/output/assembly_graph_after_simplification.gfa -c benchmark/fastq/5-strain-HIV-20000x/output/contigs.paths -mincov 500 -minlen 8000 -maxlen 10000 -overlap 127 -ref benchmark/strains/5-strain-HIV.fasta -o acc_5_hiv/ > 5hiv.log
python eval_script/quast_evaluation.py -c1 benchmark/fastq/5-strain-HIV-20000x/output/contigs.fasta -c2 acc_5_hiv/extended_contig.fasta -c3 acc_5_hiv/extended_contig.fasta -ref benchmark/strains/5-strain-HIV.fasta -o quast5hiv/

# 6 POLIO
time python src/hap_construction.py -gfa benchmark/fastq/6-strain-poliovirus/output/assembly_graph_after_simplification.gfa -c benchmark/fastq/6-strain-poliovirus/output/contigs.paths -mincov 500 -minlen 6000 -maxlen 8000 -overlap 127 -ref benchmark/strains/6-strain-polio.fasta -o acc_6polio/ > 6polio.log
python eval_script/quast_evaluation.py -c1 benchmark/fastq/6-strain-poliovirus/output/contigs.fasta -c2 ../vg-flow/vgflow6polioResult/haps.final.fasta -c3 acc_6polio/extended_contig.fasta -ref benchmark/strains/6-strain-polio.fasta -o quast6polio/

# 15 ZIKV
time python src/hap_construction.py -gfa benchmark/fastq/15-strain-ZIKV-20000x/output/assembly_graph_after_simplification.gfa -c benchmark/fastq/15-strain-ZIKV-20000x/output/contigs.paths -mincov 500 -minlen 8000 -maxlen 11000 -overlap 127 -ref benchmark/strains/15-strain-ZIKV.fasta -o acc_15_zikv/ > 15zikv.log
python eval_script/quast_evaluation.py -c1 benchmark/fastq/15-strain-ZIKV-20000x/output/contigs.fasta -c2 ../vg-flow/vgflow15o/haps.final.fasta -c3 acc_15_zikv/extended_contig.fasta -ref benchmark/strains/15-strain-ZIKV.fasta -o quast15zikv/