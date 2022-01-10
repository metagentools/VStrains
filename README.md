## The project is aiming to construct full-length haplotype from metagenomic environment, using pair-end reads

# command to assemble pair-end reads using metaviral:
```sh
$ spades --metaviral -1 forward.fastq -2 reverse.fastq -o dir
```

1. for all neg-neg edge, reverse to pos-pos
2. Pick highest depth node as source, start BFS
3. within BFS, fix any bubble node (both pos-neg been used) via estimated edge depth and longest path the node can reach)