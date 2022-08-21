#!/usr/bin/python
import subprocess
import argparse
import sys
import os

usage = "Use MetaQUAST to evaluate assembly result"
Author = "Runpeng Luo"


def sep_ref(ref_file, id=0):
    ref_file_list = []
    i = 0
    with open(ref_file, "r") as ref:
        j = 0
        lines = ref.readlines()
        l = len(lines)
        while j < l - 1:
            name_in_file = lines[j]
            name = str(lines[j][1:-1])
            name = name.split(" ")[0]
            name = name.split(".")[0]
            strain = lines[j + 1]
            j = j + 2
            file_name = "sub_" + str(id) + "_" + str(name) + "_ref.fasta"
            subprocess.check_call("touch {0}".format(file_name), shell=True)
            with open(file_name, "w") as sub_file:
                sub_file.write(name_in_file)
                sub_file.write(strain)
                sub_file.close()
            ref_file_list.append(file_name)
            i = i + 1
        ref.close()
    print("ref list: ", ref_file_list)
    return ref_file_list


def quast_eval(files, ref, o, quast, id=0):
    subprocess.check_call("rm -rf sub_{0}_*_ref.fasta".format(id), shell=True)

    ref_file_list = sep_ref(ref, id)

    command = "python2 {0} --unique-mapping -m 250 -t 8 ".format(quast)

    for fname in files:
        command += fname + " "

    command += "-o " + o + " -R "

    for file in ref_file_list:
        command += file + ","
    command = command[:-1]

    print(command)
    subprocess.check_call(command, shell=True)

    # clean up
    subprocess.check_call("rm -rf sub_{0}_*_ref.fasta".format(id), shell=True)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="quast_evaluation.py", description=usage)
    parser.add_argument(
        "-quast",
        "--path_to_quast",
        dest="quast",
        required=True,
        help="path to MetaQuast python script",
    )
    parser.add_argument(
        "-cs",
        "--contig_files",
        dest="files",
        default=None,
        nargs="+",
        help="contig files from different tools, separated by space",
    )
    parser.add_argument(
        "-d",
        "--contig_dir",
        dest="idir",
        default=None,
        help="contig files from different tools, stored in the directory",
    )
    parser.add_argument(
        "-ref",
        "--ref_file",
        dest="ref_file",
        type=str,
        required=True,
        help="ref file (single)",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        type=str,
        required=True,
        help="output directory",
    )
    args = parser.parse_args()

    if args.idir != None and (
        not os.path.exists(args.idir) or not os.path.isdir(args.idir)
    ):
        print("Please provide correct directory")
        sys.exit(1)

    if (args.idir == None and args.files == None) or (
        args.idir != None and args.files != None
    ):
        print("Please provide correct query input")
        sys.exit(1)
    files = (
        args.files
        if args.files != None
        else [str(args.idir) + s for s in sorted(os.listdir(args.idir))]
    )

    quast_eval(files, args.ref_file, args.output_dir, args.quast)
