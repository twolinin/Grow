
SD=./StriDe/grow

thread=$1
TGSPrefix=$2
ContigPrefix=$3


$SD index -t $thread $TGSPrefix.fasta
$SD index -t $thread $ContigPrefix.fasta


$SD reassembler --first  -t $thread -p $TGSPrefix -P $ContigPrefix $ContigPrefix.fasta
$SD reassembler --second -t $thread -p $TGSPrefix -P $ContigPrefix Reassembler_first_Overlap.fa
$SD reassembler --third  -t $thread -p $TGSPrefix -P $ContigPrefix Reassembler_second_NonOverlap.fa


