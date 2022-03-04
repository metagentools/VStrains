#!/usr/bin/env python3
import subprocess

if __name__ == "__main__":
    # remove all the acc directory
    subprocess.check_call("rm -rf eval_result && mkdir eval_result/", shell=True)

    # 5 strain hiv 20000
    comm1 = """time python src/hap_construction.py -gfa benchmark/fastq/5-strain-HIV-20000x/output/assembly_graph_after_simplification.gfa \
    -c benchmark/fastq/5-strain-HIV-20000x/output/contigs.paths -mincov 500 -minlen 8000 -maxlen 10000 -overlap 127 \
    -ref benchmark/strains/5-strain-HIV.fasta -o acc_5hiv_20000/ > hap_5hiv_20000.log
    """
    subprocess.check_call(comm1, shell=True)

    # 6 strain polio
    comm2 = """time python src/hap_construction.py -gfa benchmark/fastq/6-strain-poliovirus/output/assembly_graph_after_simplification.gfa \
    -c benchmark/fastq/6-strain-poliovirus/output/contigs.paths -mincov 500 -minlen 6000 -maxlen 8000 -overlap 127 \
    -ref benchmark/strains/6-strain-polio.fasta -o acc_6polio/ > hap_6polio.log
    """
    subprocess.check_call(comm2, shell=True)

    # 10 strain hcv 20000
    comm3 = """time python src/hap_construction.py -gfa benchmark/fastq/10-strain-HCV-20000x/output/assembly_graph_after_simplification.gfa \
    -c benchmark/fastq/10-strain-HCV-20000x/output/contigs.paths -mincov 500 -minlen 7000 -maxlen 9500 -overlap 127 \
    -ref benchmark/strains/10-strain-HCV.fasta -o acc_10hcv_20000/ > hap_10hcv_20000.log
    """
    subprocess.check_call(comm3, shell=True)

    # 15 strain zikv 20000
    comm4 = """time python src/hap_construction.py -gfa benchmark/fastq/15-strain-ZIKV-20000x/output/assembly_graph_after_simplification.gfa \
    -c benchmark/fastq/15-strain-ZIKV-20000x/output/contigs.paths -mincov 500 -minlen 8000 -maxlen 11000 -overlap 127 \
    -ref benchmark/strains/15-strain-ZIKV.fasta -o acc_15zikv_20000/ > hap_15zikv_20000.log
    """
    subprocess.check_call(comm4, shell=True)

    # 15 strain zikv 20000 careful
    comm5 = """time python src/hap_construction.py -gfa benchmark/fastq/15-strain-ZIKV-20000x/output_careful/assembly_graph_after_simplification.gfa \
    -c benchmark/fastq/15-strain-ZIKV-20000x/output_careful/contigs.paths -mincov 500 -minlen 8000 -maxlen 11000 -overlap 127 \
    -ref benchmark/strains/15-strain-ZIKV.fasta -o acc_15zikv_20000_careful/ > hap_15zikv_20000_careful.log
    """
    subprocess.check_call(comm5, shell=True)

    # relocate all the log file
    subprocess.check_call("mv hap_*.log eval_result/ && mv acc_* eval_result/", shell=True)