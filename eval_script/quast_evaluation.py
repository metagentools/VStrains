#!/usr/bin/python
import subprocess
import argparse

usage = "Use meta-Quast to evaluate assembly result"
Author = "Runpeng Luo"

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
    print("ref list: ", ref_file_list)
    return ref_file_list

def quast_eval(c1=None, c2=None, c3=None, ref=None, o=None, id=0):
    assert c1 and c2 and c3, ref and o
    subprocess.check_call("rm -rf sub_{0}_*_ref.fasta".format(id), shell=True)

    ref_file_list = sep_ref(ref, id)

    command = "python2 /Users/luorunpeng/bio_tools/quast-5.1.0rc1/metaquast.py --unique-mapping -m 250 {0} {1} {2} -t 8 -o {3} -R ".format(c1, c2, c3, o)

    for file in ref_file_list:
        command = command + file + ","
    command = command[:-1]
    print(command)
    subprocess.check_call(command, shell=True)

    # clean up
    subprocess.check_call("rm -rf sub_{0}_*_ref.fasta".format(id), shell=True)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='quast_evaluation.py', description=usage)
    parser.add_argument('-c1', '--contig_file1', dest='contig_file1', type=str, required=True, help='contig file 01')
    parser.add_argument('-c2', '--contig_file2', dest='contig_file2', type=str, required=True, help='contig file 02')
    parser.add_argument('-c3', '--contig_file3', dest='contig_file3', type=str, required=True, help='contig file 02')
    parser.add_argument('-ref', '--ref_file', dest='ref_file', type=str, required=True, help='ref file')
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str, required=True, help='output_dir')
    args = parser.parse_args()
    assert args.contig_file1 and args.contig_file2 and args.contig_file3 and args.ref_file and args.output_dir

    quast_eval(args.contig_file1, args.contig_file2, args.contig_file3, args.ref_file, args.output_dir)


