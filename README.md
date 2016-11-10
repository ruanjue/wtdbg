# wtdbg

A fuzzy Bruijn graph (FBG) approach to long noisy reads assembly

**wtdbg** is desiged to assemble huge genomes in very limited time, it requires a PowerPC with multiple-cores and very big RAM (1Tb+).
**wtdbg** can assemble a 100 X human pacbio dataset within one day.

# Assembly of big genomes using PacBio data

```sh
# assembly of contigs
wtdbg-1.1.006 -t 96 -i pb-reads.fa -o dbg -H -k 21 -S 1.02 -e 3 2>&1 | tee log.wtdbg
# -t: number of threads, please type 'wtdbg-1.1.006 -h' to get a document
# -i: you can set more than one sequences files, such as -i 1.fa. -i 2.fq -i 3.fa.gz -i 4.fq.gz
# -o: the prefix of results
# -S: 1.01 will use all kmers, 1.02 will use half by sumsampling, 1.04 will use 1/4, and so on
#     2.01 will use half by picking minimizers, but not fully tested
# -e: if too low coverage(< 30 X), try to set -e 2
# please note that dbg.ctg.fa is full of errors from raw reads

# first round of polishment
wtdbg-cns -t 96 -i dbg.ctg.lay -o dbg.ctg.lay.fa -k 15 2>&1 | tee log.cns.1
# dbg.ctg.lay.fa is the polished contigs

# if possible, further polishment
minimap -t 96 -L 100 dbg.ctg.lay.fa pb-reads.fa 2> >(tee log.minimap) | best_minimap_hit.pl | awk '{print $6"\t"$8"\t"$9"\t"$1"\t"$5"\t"$3"\t"$4}' >dbg.map
map2dbgcns dbg.ctg.lay.fa pb-reads.fa dbg.map >dbg.map.lay
wtdbg-cns -t 96 -i dbg.map.lay -o dbg.map.lay.fa -k 13 2>&1 | tee log.cns.2
# you need to concat all reads into one file for minimap and map2dbgcns
# dbg.map.lay.fa is the final contigs

```

# Assembly of huge genomes

**wtdbg** can handle huge genomes (10G+). Test showed that, for 20X PacBio data of a plant genome of 10G bp, contigs can be assembled within one day.


# Contact
Jue Ruan <ruanjue@gmail.com>
