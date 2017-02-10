#!/usr/bin/perl -w

$species=$ARGV[0];
$analysis=$ARGV[1];
$base_dir=$ARGV[2];

open LOG, ">>${base_dir}/Analysis/log.txt";

open OUTPUT, ">${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_rbs_coordinates.tab";

open INPUT, "${base_dir}/Data/Reference_files/${species}_reference.fasta";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^([ATGCN]+)/){
		$seq=$1;
		@seq_array=split(//, $seq);
	}
}

open INPUT, "${base_dir}/Analysis/${analysis}/${species}_${analysis}/${species}_rbs_finder_output.tab";
while(<INPUT>){
	$line=$_;
	chomp $line;
	
	if($line =~ /^\s*(\S+\s+\d+\s+\d+.*)/){
		$line=$1;
		@line_array=split(/\s+/, $line);
		
		$gene_id=$line_array[0];
		$new_sta=$line_array[1];
		$old_sta=$line_array[8];
		$pat=$line_array[3];
		$pos=$line_array[4];
		
		$ind_pos=$pos-1;
		$pat_len=length($pat);
		
		$seq_pat="";
		if($new_sta == $old_sta && $pat_len == 5){
			if($pos < $old_sta){
				$end=$pos+4;
				
				for($i=$ind_pos; $i<($ind_pos+5); $i++){
					$seq_pat="$seq_pat$seq_array[$i]";
				}
				
				if($pat eq $seq_pat){
					print OUTPUT "$gene_id\t$pos\t$end\n";
				}else{
					print "bad pattern\n";
				}
			}elsif($pos > $old_sta){
				$pos=$ind_pos-4;
				$ind_pos=$ind_pos-4;
				$end=$pos+4;
				
				for($i=$ind_pos; $i<($ind_pos+5); $i++){
					$seq_pat="$seq_pat$seq_array[$i]";
				}
				
				$seq_pat=~tr/ATGC/TACG/;
				$seq_pat=reverse($seq_pat);
				
				if($pat eq $seq_pat){
					print OUTPUT "$gene_id\t$pos\t$end\n";
				}else{
					print "bad pattern\n";
				}
			}
		}
	}
}

print "$species rbs predictions completed.\n";

print LOG "$species rbs predictions completed.\n";

