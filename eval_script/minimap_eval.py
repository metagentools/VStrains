import subprocess
import argparse
def map_ref_to_graph(paf_file):
    """
    map reference strain to the graph, debug only
    assumption: graph is stored in acc/simplifed_graph, 
    """
    strain_dict = {}
    all = set()
    with open(paf_file, 'r') as paf:
        for Line in paf:
            splited = Line.split('\t')
            seg_no = str(splited[0])
            seg_no = splited[0]
            seg_l = int(splited[1])
            seg_s = int(splited[2])
            seg_f = int(splited[3])
            ref_no = str(splited[5])
            nmatch = int(splited[9])
            nblock = int(splited[10])
            mark = int(splited[11])
            if ref_no not in strain_dict:
                strain_dict[ref_no] = []
            all.add(seg_no)
            strain_dict[ref_no].append((seg_no,nmatch/nblock))
        paf.close()
    
    print("strain dict mapping")
    max_mappings = set()
    for ref_no, cands in strain_dict.items():
        sorted_cands = sorted(cands, key=lambda x:x[1], reverse=True)
        print("Strain: ", ref_no, " Cands: ", sorted_cands)
        for candno, candr in sorted_cands:
            if candno not in max_mappings:  
                print(ref_no, candno)
                max_mappings.add(candno)
                break
        print("-------------------")
    print("Size: ", len(max_mappings), "max mappings: ", max_mappings)
    print("DUP: ", all.difference(max_mappings))
    return strain_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='minimap_eval.py')
    parser.add_argument('-paf', '--minimap_file', dest='paf_file', type=str, required=True, help='minimap paf file')
    args = parser.parse_args()
    map_ref_to_graph(args.paf_file)
