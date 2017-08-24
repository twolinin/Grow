# Introduction
The StriDe Assembler integrates string and de Bruijn graph by decomposing reads within error-prone regions, while extending paire-end read into long reads for assembly through repetitive regions. The entire implementation is done by revising Simpson's [SGA][1] components for our own purpose, porting Li's [ropebwt2][2] for FM-index construction, and adding key components of this assembler. 

# Executable version
A precompiled version under Linux 64bit (StriDe_Linux_64bit) can be directly downloaded and executed. 

# Compile by yourself
To compile StriDe assembler in your specific environment, type 

      1. ./autogen.sh 
      2. ./configure
      3. make

An executable program called stride will be found under the StriDe folder.

# Execution

	stride index -t 10 ContigPrefix.fasta
	stride index -t 10 SeqReadPrefix.fasta

	stride reassembler --first  -t 10 -p SeqReadPrefix -P ContigPrefix ContigPrefix.fasta
	stride reassembler --second -t 10 -p SeqReadPrefix -P ContigPrefix Reassembler_first_Overlap.fa
	stride reassembler --third  -t 10 -p SeqReadPrefix -P ContigPrefix Reassembler_second_NonOverlap.fa
