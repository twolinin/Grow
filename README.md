# Grow
Grow: Genome Reassembly at Low-Coverage Regions for Third-Generation Sequencing

Third-generation sequencing (TGS) is becoming the preferred choice for de novo genome assembly mainly due to much longer reads and fast turnaround time. Unfortunately, existing TGS assemblers often failed to assemble through low-coverage regions because of insufficient sequencing depth or uneven coverage across the entire genome. This is further complicated by the overlapping algorithm underlying each assembler, which often employ dimensional reduction for speedup but sacrifice sensitivity in these regions.

# Compile by yourself
To compile StriDe assembler in your specific environment, type 

      1. ./autogen.sh 
      2. ./configure
      3. make

An executable program called grow will be found under the StriDe folder.

# Execution

Example of a Step-by-step script. The Contigs are assembled by raw reads, and the raw read must be third generation sequencing data.

	grow index -t 10 Contig.fasta
	grow index -t 10 RawRead.fasta

	grow reassembler --first  -t 10 -p RawRead -P Contig Contig.fasta
	grow reassembler --second -t 10 -p RawRead -P Contig Reassembler_first_Overlap.fa
	grow reassembler --third  -t 10 -p RawRead -P Contig Reassembler_second_NonOverlap.fa


# Related Links
https://github.com/ythuang0522/StriDe

