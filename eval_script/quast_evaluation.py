#!/usr/bin/python
from re import sub
import subprocess
import argparse

usage = "separate ref"
def sep_ref(ref_file):
    ref_file_list = []
    i = 0
    with open(ref_file, 'r') as ref:
        j = 0
        lines = ref.readlines()
        l = len(lines)
        while j < l - 1:
            name = lines[j]
            strain = lines[j+1]
            j = j + 2
            file_name = "sub_" + str(name[1:-1]) + "_ref.fasta"
            subprocess.check_call("echo "" > {0}".format(file_name), shell=True)
            with open (file_name, 'w') as sub_file:
                sub_file.write(name)
                sub_file.write(strain)
                sub_file.close()
            ref_file_list.append(file_name)
            i = i + 1
        ref.close()
    return ref_file_list

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='quast_query.py', description=usage)
    parser.add_argument('-c', '--contig_file', dest='contig_file', type=str, required=True, help='contig file')
    parser.add_argument('-ref', '--ref_file', dest='ref_file', type=str, required=True, help='ref file')
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str, required=True, help='output_dir')
    args = parser.parse_args()

    subprocess.check_call("rm -rf sub_*_ref.fasta", shell=True)

    ref_file_list = sep_ref(args.ref_file)
    
    command = "/Users/luorunpeng/bio_tools/quast-5.1.0rc1/metaquast.py {0} -o {1} -R ".format(args.contig_file, args.output_dir)
    for file in ref_file_list:
        command = command + file + ","
    command = command[:-1]
    subprocess.check_call(command, shell=True)

    # clean up
    subprocess.check_call("rm -rf sub_*_ref.fasta", shell=True)
