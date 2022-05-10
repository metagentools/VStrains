
paths = [[1, 3, 4, 3, 1, 2], 
         [1, 3, 3, 3, 2, 2]]
contig_list = []
k = -1
for i in range(6):
    currnodes = [path[i] for path in paths]
    if all(e==currnodes[0] for e in currnodes):
        if k == i - 1:
            # just concat to the last sublist from contig_list
            if len(contig_list) == 0:
                contig_list = [[currnodes[0]]]
            else:
                contig_list[-1].append(currnodes[0])
            k += 1
        else:
            contig_list.append([currnodes[0]])
            k = i
print(contig_list)