#!/usr/bin/python
from re import sub
import subprocess
import argparse

usage = "separate ref"
def sep_ref(ref_file, id=0):
    ref_file_list = []
    i = 0
    with open(ref_file, 'r') as ref:
        j = 0
        lines = ref.readlines()
        l = len(lines)
        while j < l - 1:
            name_in_file = lines[j]
            name = str(lines[j][1:-1])
            name = name.split(" ")[0]
            name = name.split(".")[0]
            strain = lines[j+1]
            j = j + 2
            file_name = "sub_" + str(id) + "_" + str(name) + "_ref.fasta"
            subprocess.check_call("touch {0}".format(file_name), shell=True)
            with open (file_name, 'w') as sub_file:
                sub_file.write(name_in_file)
                sub_file.write(strain)
                sub_file.close()
            ref_file_list.append(file_name)
            i = i + 1
        ref.close()
    return ref_file_list

def quast_eval(c=None, ref=None, o=None, id=0):
    assert c and ref and o
    subprocess.check_call("rm -rf sub_{0}_*_ref.fasta".format(id), shell=True)

    ref_file_list = sep_ref(ref, id)

    command = "/Users/luorunpeng/bio_tools/quast-5.1.0rc1/metaquast.py -m 100 {0} -t 8 -o {1} -R ".format(c, o)

    for file in ref_file_list:
        command = command + file + ","
    command = command[:-1]
    subprocess.check_call(command, shell=True)

    # clean up
    subprocess.check_call("rm -rf sub_{0}_*_ref.fasta".format(id), shell=True)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='quast_query.py', description=usage)
    parser.add_argument('-c', '--contig_file', dest='contig_file', type=str, required=True, help='contig file')
    parser.add_argument('-ref', '--ref_file', dest='ref_file', type=str, required=True, help='ref file')
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str, required=True, help='output_dir')
    args = parser.parse_args()
    assert args.contig_file and args.ref_file and args.output_dir

    quast_eval(args.contig_file, args.ref_file, args.output_dir)


