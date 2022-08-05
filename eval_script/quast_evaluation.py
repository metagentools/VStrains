#!/usr/bin/python
import subprocess
import argparse

usage = "Use meta-Quast to evaluate assembly result"
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


def quast_eval(files, ref=None, o=None, id=0):
    subprocess.check_call("rm -rf sub_{0}_*_ref.fasta".format(id), shell=True)

    ref_file_list = sep_ref(ref, id)

    command = "python2 /Users/luorunpeng/bio_tools/quast-5.1.0rc1/metaquast.py --unique-mapping -m 250 -t 8 "

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
        "-cs",
        "--contig_files",
        dest="files",
        nargs="+",
        required=True,
        help="<Required> contig files",
    )
    parser.add_argument(
        "-ref", "--ref_file", dest="ref_file", type=str, required=True, help="ref file"
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        type=str,
        required=True,
        help="output_dir",
    )
    args = parser.parse_args()

    quast_eval(args.files, args.ref_file, args.output_dir)
