#!/usr/bin/ perl
#############################################################################
#
# by Baris Ethem Suzek suzek@cs.jhu.edu 
#    12/19/2001
#    rbs_finder.pl 
#
# Usage: 
#      perl rbs_finder.pl <whole_dna_sequence> <whole_gene_coordinate_list> <output_file> <rbs_region_length> [<consensus_sequence>] [-p <to_be_relocated_gene_coord_list>] 
#      
#      where:
#           whole_dna_sequence: the file containning all the sequence
#           to_be_relocated_gene_coordinate_list: the file in the format "GeneId StartCoord EndCoord" and 
#                                                 which contains the list of gene coords to be relocated
#           whole_gene_coordinate_list: the file in the format "GeneId StartCoord EndCoord" and contains
#                                       all the gene coords 
#           output_file: the name of the output file
#           rbs_region_length: the window size to look RBS for
#           consensus_sequence:(optional) the consensus sequence to look for in the upstream regions of start codons
# 
# NOTE: A low scoring ATG is prefered over a high scoring GTG, a low scoring ATG or GTG is prefered 
# over a high scoring TTG.
#
#############################################################################

open SEQ_FILE,$ARGV[0];
open OUTPUT , ">$ARGV[2]";
$RBS_region_len=$ARGV[3];
#initilaize consensus sequence

if (($ARGV[4] ne "") && ($ARGV[4] !~ /\-p/)){ 
    $SD_seq=lc($ARGV[4]);
    if ($SD_seq !~ /[atgc]+/){
	print "Error in given consensus sequence\n";
	exit;
    }
    if ($ARGV[5] =~ /\-p/){
	if ($ARGV[6] eq "") {print "Error in partial gene list file\n"; exit;}
	$gene_coord_file=$ARGV[6];
    }
    else{
	$gene_coord_file=$ARGV[1];
    }
} elsif (($ARGV[4] =~ /\-p/)){ 
    $SD_seq="aggag";
    if ($ARGV[5] eq "") {print "Error in partial gene list file\n"; exit;}
    $gene_coord_file=$ARGV[5];

} else {
    $SD_seq="aggag";
    $gene_coord_file=$ARGV[1];
}


#read the whole sequence
$seq="";
while (<SEQ_FILE>){
	chomp;
	$_ = lc($_);
	s/~[a-z]//g;
	next if ($_=~/^>/);
	$seq.=lc($_);
}

&init_operons();
&train_SD_markov();



print OUTPUT "                NewStart     Stop      RibosomeBindingSite                                OldStart \n";
print OUTPUT "    GeneID      Position    Position     Pattern Position     NewStartCodon    Shift   Codon Position\n"; 

$rbs_relocated=0;
$rbs_original=0;
$no_rbs=0;

#****************************
open COORD_FILE,$gene_coord_file;
while(<COORD_FILE>){
	chomp;
	next if ($_=~ /^\#/);
	next if ($_  eq ""); 
	
        #**** glimmer2 before get-putative
	if ($_ =~ /GC Proportion/){
	    while(1){
		$skip=<COORD_FILE>;
		last if ($skip=~ /Putative\sGenes/);
	    }
	    next;
	}
	#*********************************

	@line_arr=split(/\s+/,$_);
	#for glimmer format
	if ($line_arr[0] eq ""){
		$start_coord=$line_arr[2];
		$end_coord=$line_arr[3];
		$orf_name=$line_arr[1];
	}
	#for <geneid> <startcoord> <endcoord> 
	else{
		$start_coord=$line_arr[1];
		$end_coord=$line_arr[2];
		$orf_name=$line_arr[0];
	}
#	print "start = $start_coord, stop = $end_coord\n";
	$RBS="";
	$max_RBS_score=-100000000;
	$max_RBS_loc=-1;
	$max_RBS_start=-1;
	$max_RBS="";
	if ($start_coord<$end_coord) {
	    $original_start_codon=substr($seq,$start_coord-1,3);
	    $new_start_upstream=&next_start_coord_upstream($original_start_codon);
	    ($new_start_coord_up,$RBS_up,$RBS_loc_up,$RBS_score_up)=split(" ",$new_start_upstream);
#	    if ($new_start_upstream eq ""){
	    $how_far_to_move=int(($end_coord-$start_coord)*(($end_coord-$start_coord+1)>300?0.35:0.2));
	    if ($how_far_to_move>300) { $how_far_to_move=300};
	    $how_far_to_move=$start_coord+$how_far_to_move;
	    $new_start_coord=$start_coord-1;
	    while($new_start_coord<$how_far_to_move) {
		if ($new_start_coord>$RBS_region_len){ $RBS_region=substr($seq,$new_start_coord-$RBS_region_len,$RBS_region_len);}
		else{ $RBS_region=substr($seq,0,$new_start_coord);}
		$tmp_threshold=&calculate_threshold(abs($how_far_to_move-$start_coord),($new_start_coord-$start_coord+1));	
		($RBS,$RBS_loc,$RBS_score)=split(" ",&check_RBS_SD($RBS_region,$tmp_threshold));
		$codon1=substr($seq,$new_start_coord,3);
		$codon2=substr($seq,$max_RBS_start,3);
		if (($RBS ne "") && ((&sc_start($codon1)>&sc_start($codon2))
		    || (($RBS_score>$max_RBS_score) 
			&& (&sc_start($codon1)==&sc_start($codon2))))){
		    $max_RBS=$RBS;
		    $max_RBS_score=$RBS_score;
		    $max_RBS_loc=$RBS_loc;
		    $max_RBS_start=$new_start_coord;
#		    last if ($new_start_coord==($start_coord-1) && ($RBS_score>$initial_threshold)) ;
		}
		$new_start_coord=&next_start_coord_fw($new_start_coord,$original_start_codon);
	    }
	    
	    $codon1=substr($seq,$new_start_coord_up,3);
	    $codon2=substr($seq,$max_RBS_start,3);
	    if  (($new_start_upstream  ne "" )
		 && ((($RBS_score_up>$max_RBS_score) 
		      && (&sc_start($codon1)==&sc_start($codon2)))
		     || (&sc_start($codon1)>&sc_start($codon2)))){
		printf(OUTPUT "%10s      %8d    %8d    %8s %8d        %5s         %5d   %5s %8d\n",
		       $orf_name,
		       $new_start_coord_up+1,
		       $end_coord,
		       uc($RBS_up),
		       $new_start_coord_up-$RBS_region_len+$RBS_loc_up+1,
		       uc(substr($seq,$new_start_coord_up,3)),
		       ($new_start_coord_up+1-$start_coord),
		       uc($original_start_codon),
		       $start_coord);
		$rbs_relocated++;
		$upstreamreloc.=" ".$orf_name;
		$upstream++;
	    }
	    else{
		if ($max_RBS ne ""){
		    printf(OUTPUT "%10s      %8d    %8d    %8s %8d        %5s         %5d   %5s %8d\n",
			   $orf_name,
			   $max_RBS_start+1,
			   $end_coord,
			   uc($max_RBS),
			   $max_RBS_start-$RBS_region_len+$max_RBS_loc+1,
			   uc(substr($seq,$max_RBS_start,3)),
			   abs($max_RBS_start+1-$start_coord),
			   uc($original_start_codon),
			   $start_coord);
		    ($start_coord == ($max_RBS_start+1))?($rbs_original++):($rbs_relocated++);
		}
		else{
		    printf(OUTPUT "%10s      %8d    %8d    %8s %8d        %5s         %5d   %5s %8d\n",
			   $orf_name,
			   $start_coord,
			   $end_coord,
			   "---",
			   "---",
			   uc($original_start_codon),
			   0,
			   uc($original_start_codon),
			   $start_coord);
		    $no_rbs++;
		}
	    }
	}   
	else{ #if the orf is on the reverse strain
	    $original_start_codon=&reverse_complement(substr($seq,$start_coord-3,3));
	    $new_start_upstream=&next_start_coord_upstream($original_start_codon);
	    ($new_start_coord_up,$RBS_up,$RBS_loc_up,$RBS_score_up)=split(" ",$new_start_upstream);
	    $how_far_to_move=int(($start_coord-$end_coord)*(abs($end_coord-$start_coord+1)>300?0.35:0.2));
	    if ($how_far_to_move>300) {$how_far_to_move=300;}
	    $how_far_to_move=$start_coord-$how_far_to_move;
	    $new_start_coord=$start_coord-1;
	    while($new_start_coord>$how_far_to_move){
		if (($new_start_coord+$RBS_region_len)<length($seq)){$RBS_region=&reverse_complement(substr($seq,$new_start_coord+1,$RBS_region_len));}
		else{$RBS_region=&reverse_complement(substr($seq,$new_start_coord+1));}
		$tmp_threshold=&calculate_threshold(abs($how_far_to_move-$start_coord),($start_coord-$new_start_coord-1)); 
		($RBS,$RBS_loc,$RBS_score)=split(" ",&check_RBS_SD($RBS_region,$tmp_threshold));
		$codon1=substr($seq,$new_start_coord-2,3);
		$codon2=substr($seq,$max_RBS_start-2,3);
		if(($RBS ne "") 
		   && ((&sc_start($codon1)>&sc_start($codon2))
		       || (($RBS_score>$max_RBS_score) 
			   && (&sc_start($codon1) == &sc_start($codon2))))){
		    $max_RBS=$RBS;
		    $max_RBS_score=$RBS_score;
		    $max_RBS_loc=$RBS_loc;
		    $max_RBS_start=$new_start_coord;
#		    last if ($new_start_coord==($start_coord-1) && ($RBS_score>$initial_threshold)) ;
		}		    
		$new_start_coord=&next_start_coord_reverse($new_start_coord,$original_start_codon);
	    }
	    $codon1=substr($seq,$new_start_coord_up-2,3);
            $codon2=substr($seq,$max_RBS_start-2,3);
	    if (($new_start_upstream ne "") 
		&& ((($RBS_score_up>$max_RBS_score) 
		     && (&sc_start($codon1)==&sc_start($codon2))) 
		    || (&sc_start($codon1)>&sc_start($codon2)))){
		    printf(OUTPUT "%10s      %8d    %8d    %8s %8d        %5s         %5d   %5s %8d\n",
		       $orf_name,
		       $new_start_coord_up+1,
		       $end_coord,
		       uc($RBS_up),
		       $new_start_coord_up+1+$RBS_region_len-$RBS_loc_up,   
		       uc(&reverse_complement(substr($seq,$new_start_coord_up-2,3))),
		       -1*($new_start_coord_up+1-$start_coord),
		       uc($original_start_codon),
		       $start_coord,);
		$rbs_relocated++;
		$upstream++;
		$upstreamreloc.=" ".$orf_name;
	    }
	    else{
		if ($max_RBS ne ""){
		    printf(OUTPUT "%10s      %8d    %8d    %8s %8d        %5s         %5d   %5s %8d\n",
			   $orf_name,
			   $max_RBS_start+1,
			   $end_coord,
			   uc($max_RBS),
			   $max_RBS_start+1+$RBS_region_len-$max_RBS_loc,   
			   uc(&reverse_complement(substr($seq,$max_RBS_start-2,3))),
			   abs($max_RBS_start+1-$start_coord),
			   uc($original_start_codon),
			   $start_coord);
		    ($start_coord == ($max_RBS_start+1))?($rbs_original++):($rbs_relocated++);
		}
		else{
		    printf(OUTPUT "%10s      %8d    %8d    %8s %8d        %5s         %5d   %5s %8d\n",
			   $orf_name,
			   $start_coord,
			   $end_coord,
			   "---",
			   "---",
			   uc($original_start_codon),
			   0,   
			   uc($original_start_codon),
			   $start_coord);
		    $no_rbs++;
		}
	    }			
	}		
    }
    

$total_orfs=$no_rbs+$rbs_original+$rbs_relocated;
print "Summary:\n";
printf("# of orfs that have RBS before original start codon loc= %d \-\> %4.2f\%\n",$rbs_original,$rbs_original*100/$total_orfs); 
printf("# of orfs that have RBS before new start codon loc= %d \-\> %4.2f\%\n",$rbs_relocated,$rbs_relocated*100/$total_orfs); 
printf("# of orfs that have no RBS= %d \-\>  %4.2f\%\n",$no_rbs,$no_rbs*100/$total_orfs); 
print "Total # of orfs: $total_orfs\n";
#print "Upstream $upstream $upstreamreloc\n";

sub next_start_coord_fw
{
   my $tmp_start=shift;
    my $o_start=shift;
    while(1){
	$tmp_start+=3;
	if ($tmp_start >= length($seq)) {return length($seq)};
	my $start_codon=substr($seq,$tmp_start,3);
	if (($start_codon eq "atg") || (($start_codon eq "gtg") || ($start_codon eq "ttg")))  {
	    return $tmp_start;
	}
    }
}

sub next_start_coord_reverse
{
    my $tmp_start=shift;
    my $o_start=shift;
    while(1){
	$tmp_start-=3;
	if ($tmp_start <=0 ) {return 0};
	my $start_codon=substr($seq,$tmp_start-2,3);
	if (($start_codon eq "cat") || (($start_codon eq "cac") || ($start_codon eq "caa"))) {
	    return $tmp_start;
	}
    }
}

sub next_start_coord_upstream
{
    my $o_start=shift;
    my $tmp_RBS_region;
    my $tmp_RBS;
    my $tmp_RBS_loc;
    my $tmp_threshold=0;	
    my $codon;
    my $y;
    my $tmp_max_RBS_score=-100000000;
    my @tmp_RBS_arr=();
    my $tmp_max_RBS="";
    my $tmp_max_loc;

    if ($start_coord<$end_coord){ #the orf processed is in forward strain
	if ($orf_previous{$orf_name}>0){ #if previous orf is in same direction
	    if ($orf_previous{$orf_name}>=$start_coord){ return "";}
	    for ($y=$start_coord-4;$y>$orf_previous{$orf_name};$y-=3){
	        $codon=substr($seq,$y,3);
		if (($codon eq "tga") || (($codon eq "taa") || ($codon eq "tag"))){last;};
		if (($codon eq "atg") || ($codon eq "gtg") || ($codon eq "ttg")){
		    $tmp_RBS_region=substr($seq,$y-$RBS_region_len,$RBS_region_len);
		    $tmp_threshold= $markov_score_threshold;
		    $tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
		    if ($tmp_RBS ne ""){ 
			@tmp_RBS_arr=split(" ",$tmp_RBS);
			if ((&sc_start($codon)>&sc_start(substr($seq,$tmp_max_loc,3))) 
			    || (($tmp_RBS_arr[2]>$tmp_max_RBS_score)
				&& (&sc_start($codon) == &sc_start(substr($seq,$tmp_max_loc,3))))){
			    $tmp_max_RBS=$tmp_RBS;
			    $tmp_max_RBS_score=$tmp_RBS_arr[2];
			    $tmp_max_loc=$y;
			}
		    }
		}
	    }
 	    if ($tmp_max_RBS ne "") {return ($tmp_max_loc." ".$tmp_max_RBS." ".$tmp_max_RBS_score);}
	    return "";
	}
	else{ #if previous orf is in opposite direction
	    if ((-1*$orf_previous{$orf_name})>=$start_coord){ return "";}
	    $tmp_RBS_region=&reverse_complement(substr($seq,-1*$orf_previous{$orf_name},$RBS_region_len));
	    $tmp_threshold=$markov_score_threshold;
	    $tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
	    if ($tmp_RBS ne ""){ #if there is an RBS before the codon
		for ($y=$start_coord-4;$y>-1*$orf_previous{$orf_name};$y-=3){
		    $codon=substr($seq,$y,3);
		    if (($codon eq "tga") || (($codon eq "taa") || ($codon eq "tag"))){ last;}
		    if (($codon eq "atg") || ($codon eq "gtg") || ($codon eq "ttg")){
			$tmp_RBS_region=substr($seq,$y-$RBS_region_len,$RBS_region_len);
			$tmp_threshold=$markov_score_threshold;
			$tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
			if ($tmp_RBS ne ""){
			    @tmp_RBS_arr=split(" ",$tmp_RBS);
			    if ((&sc_start($codon)>&sc_start(substr($seq,$tmp_max_loc,3))) 
				|| (($tmp_RBS_arr[2]>$tmp_max_RBS_score)
				    && (&sc_start($codon)==&sc_start(substr($seq,$tmp_max_loc,3))))){
				$tmp_max_RBS=$tmp_RBS;
				$tmp_max_RBS_score=$tmp_RBS_arr[2];
				$tmp_max_loc=$y;
			    }
			}
		    }
		}
		if ($tmp_max_RBS ne "") {return ($tmp_max_loc." ".$tmp_max_RBS." ".$tmp_max_RBS_score);}
		return "";	
	    }
	    else{ #if there is no RBS before previous orf
		for ($y=$start_coord-4;$y>0;$y-=3){
		    $codon=substr($seq,$y,3);
		    if (($codon eq "tga") || (($codon eq "taa") || ($codon eq "tag"))){last;}
		    if (($codon eq "atg") || ($codon eq "gtg") || ($codon eq "ttg")){
			$tmp_RBS_region=substr($seq,$y-$RBS_region_len,$RBS_region_len);
			$tmp_threshold=$markov_score_threshold;
			$tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
			if ($tmp_RBS ne ""){
			    @tmp_RBS_arr=split(" ",$tmp_RBS);
			    if ((&sc_start($codon)>&sc_start(substr($seq,$tmp_max_loc,3))) 
				|| (($tmp_RBS_arr[2]>$tmp_max_RBS_score) 
				    && (&sc_start($codon)==&sc_start(substr($seq,$tmp_max_loc,3))))){
				$tmp_max_RBS=$tmp_RBS;
				$tmp_max_RBS_score=$tmp_RBS_arr[2];
				$tmp_max_loc=$y;
			    }
			}
		    }
		}
		if ($tmp_max_RBS ne "") {return ($tmp_max_loc." ".$tmp_max_RBS." ".$tmp_max_RBS_score);}
		return "";
	    }
	}
    }
    else{ #the orf processed is in reverse strain
	if ($orf_previous{$orf_name}<0){ #if previous orf is in same direction
	    if ((-1*$orf_previous{$orf_name})<=$start_coord){ return "";}
	    for ($y=$start_coord+2;$y<(-1*$orf_previous{$orf_name});$y+=3){
		$codon=&reverse_complement(substr($seq,$y-2,3));
		if (($codon eq "tga") || (($codon eq "taa") || ($codon eq "tag"))){last;}
		if (($codon eq "atg") || ($codon eq "gtg") || ($codon eq "ttg")){
		    $tmp_RBS_region=&reverse_complement(substr($seq,$y+1,$RBS_region_len));
		    $tmp_threshold=$markov_score_threshold;
		    $tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
		    if ($tmp_RBS ne ""){
			@tmp_RBS_arr=split(" ",$tmp_RBS);
			if ((&sc_start($codon)>&sc_start(&reverse_complement(substr($seq,$tmp_max_loc-2,3)))) 
			    || (($tmp_RBS_arr[2]>$tmp_max_RBS_score) 
				&& (&sc_start($codon)==&sc_start(&reverse_complement(substr($seq,$tmp_max_loc-2,3)))))){
			    $tmp_max_RBS=$tmp_RBS;
			    $tmp_max_RBS_score=$tmp_RBS_arr[2];
			    $tmp_max_loc=$y;
			}
		    }
		}
	    }
	    if ($tmp_max_RBS ne "") {return ($tmp_max_loc." ".$tmp_max_RBS." ".$tmp_max_RBS_score);}
	    return "";
	}
	else{#if previous orf is in opposite direction
		 if (($orf_previous{$orf_name})<=$start_coord){ return "";}
		 $tmp_RBS_region=substr($seq,$orf_previous{$orf_name}-$RBS_region_len-1,$RBS_region_len);
		 $tmp_threshold=$markov_score_threshold;
		 $tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
		 if ($tmp_RBS ne ""){ #if there is an RBS before the codon
		     for ($y=$start_coord+2;$y<$orf_previous{$orf_name};$y+=3){
			 $codon=&reverse_complement(substr($seq,$y-2,3));
			 if (($codon eq "tga") || (($codon eq "taa") || ($codon eq "tag"))){last;}
			 if (($codon eq "atg") || ($codon eq "gtg") || ($codon eq "ttg")){
			     $tmp_RBS_region=&reverse_complement(substr($seq,$y+1,$RBS_region_len));
			     $tmp_threshold=$markov_score_threshold;
			     $tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
			     if ($tmp_RBS ne ""){
				 @tmp_RBS_arr=split(" ",$tmp_RBS);
				 if (&sc_start($codon)>&sc_start(&reverse_complement(substr($seq,$tmp_max_loc-2,3)))
				     || (($tmp_RBS_arr[2]>$tmp_max_RBS_score)
					 && (&sc_start($codon)==&sc_start(&reverse_complement(substr($seq,$tmp_max_loc-2,3)))))){
				     $tmp_max_RBS=$tmp_RBS;
				     $tmp_max_RBS_score=$tmp_RBS_arr[2];
				     $tmp_max_loc=$y;
				 }
			     }
			 }
		     }
		     if ($tmp_max_RBS ne "") {return ($tmp_max_loc." ".$tmp_max_RBS." ".$tmp_max_RBS_score);}
		     return "";
		 }
		 else{
		     for ($y=$start_coord+2;$y<length($seq);$y+=3){
			 $codon=&reverse_complement(substr($seq,$y-2,3));
			 if (($codon eq "tga") || (($codon eq "taa") || ($codon eq "tag"))){last;}
			 if (($codon eq "atg") || ($codon eq "gtg") || ($codon eq "ttg")){
			     $tmp_RBS_region=&reverse_complement(substr($seq,$y+1,$RBS_region_len));
			     $tmp_threshold=$markov_score_threshold;
			     $tmp_RBS=&check_RBS_SD($tmp_RBS_region,$tmp_threshold);
			     if ($tmp_RBS ne ""){
				 @tmp_RBS_arr=split(" ",$tmp_RBS);
				if (&sc_start($codon)>&sc_start(&reverse_complement(substr($seq,$tmp_max_loc-2,3))) 
				    || (($tmp_RBS_arr[2]>$tmp_max_RBS_score)
					&& (&sc_start($codon)==&sc_start(&reverse_complement(substr($seq,$tmp_max_loc-2,3)))))){
				    $tmp_max_RBS=$tmp_RBS;
				    $tmp_max_RBS_score=$tmp_RBS_arr[2];
				    $tmp_max_loc=$y;
				 }
			     }
			 }
		     }
		     if ($tmp_max_RBS ne "") {return ($tmp_max_loc." ".$tmp_max_RBS." ".$tmp_max_RBS_score);}
		     return "";
		 }
	     }
	return "";
    }
}


sub reverse_complement{
    #reverse the input seq
    my $tmp_seq=shift;
    $tmp_seq=lc(reverse($tmp_seq));
    
    #take complement of reverse
    my @tmp_seq_array=split(//,$tmp_seq);
    $tmp_seq="";
    for (my $k=0;$k<=$#tmp_seq_array;$k++){
	if ($tmp_seq_array[$k] eq "a"){
		$tmp_seq.="t";
	}
	elsif ($tmp_seq_array[$k] eq "t"){
                $tmp_seq.="a";
	    }
	elsif ($tmp_seq_array[$k] eq "g"){
                $tmp_seq.="c";
        }
	elsif ($tmp_seq_array[$k] eq "c"){
                $tmp_seq.="g";
        }
	else{$tmp_seq.=$tmp_seq_array[$k];};
    }
    return $tmp_seq;
}


sub calculate_markov_score
{
    my $tmp_seq=shift;
    my $tmp_score=0;
    for (my $r=0;$r<length($tmp_seq);$r++){
	$tmp_score+=$SD_markov_score[$r]{substr($tmp_seq,$r,1)};
    }
    return $tmp_score;
}

sub score_tuple
{
    my $f=shift;
    my $s=shift;
    if ($f eq $s){
	return (($f eq "g" || $f eq "c")?3:2);
    }
    elsif (($f eq "a") && ($s eq "g")){
	return 2;
    }
    elsif (($f eq "c") && ($s eq "t")){
	return 2;
    }
    else {return 0};
}

 

sub train_SD_markov
{
    open ORF_LIST,$gene_coord_file;
    my $threshold_score=(length($SD_seq)-3)*3+3;
    my $number_sample=0;
    my @lowest_scored_seq;
    
	
    while(<ORF_LIST>){
	chomp;
	next if ($_=~ /^\#/);
	next if ($_  eq "");
	#**** glimmer2 before get-putative
	if ($_ =~ /GC Proportion/){
	    while(1){
		$skip=<ORF_LIST>;
		last if($skip=~ /Putative\sGenes/);
	    }
	    next;
	}
	#*********************************
	my @line_arr=split(/\s+/,$_);
	my $test_RBS_region="";
	if ($line_arr[0] eq ""){
	    if ($line_arr[2]>$line_arr[3]){
		$test_RBS_region=&reverse_complement(substr($seq,$line_arr[2],$RBS_region_len));
	    }
	    else{
		$test_RBS_region=substr($seq,$line_arr[2]-$RBS_region_len-1,$RBS_region_len);
	    }
	}
	else{
	    if ($line_arr[1]>$line_arr[2]){
	        $test_RBS_region=&reverse_complement(substr($seq,$line_arr[1],$RBS_region_len));
	    }
	    else{
		$test_RBS_region=substr($seq,$line_arr[1]-$RBS_region_len-1,$RBS_region_len);
	    }
	}
	
	next if (length($test_RBS_region) < length($SD_seq));

	#calculate all the 6 substrs of RBS region find the most 
	my $max_score=0;
	for (my $q=0;$q<=(length($test_RBS_region)-length($SD_seq));$q++){
	    my $tmp_score=0;
	    my $tmp_seq=substr($test_RBS_region,$q,length($SD_seq));
	    for (my $w=0;$w<length($SD_seq);$w++){
		$tmp_score+=&score_tuple(substr($SD_seq,$w,1),substr($tmp_seq,$w,1));
	    }
            if ($tmp_score>$max_score){
		$max_seq=$tmp_seq;
		$max_score=$tmp_score;
		$max_loc=$q;
	    }			       
	}

	#all the scores for subsequences are calculated get the max scored seq
	#to use in 0th order Markov model    

	if ($max_score>$threshold_score){
	    $RBS_seq[$number_of_sample]=$max_seq;
	    $number_of_sample++;
	    for (my $q=0;$q<length($max_seq);$q++){ 
		$SD_markov_score[$q]{substr($max_seq,$q,1)}++
	    };
	    $loc_hash{$max_loc}++;
	}    
    }
    $number_of_sample=($number_of_sample==0)?1:$number_of_sample;
    for (my $q=0;$q<length($SD_seq);$q++){
	$SD_markov_score[$q]{"a"}=log(($SD_markov_score[$q]{"a"}+0.000001)/$number_of_sample);	
	$SD_markov_score[$q]{"t"}=log(($SD_markov_score[$q]{"t"}+0.000001)/$number_of_sample);	
	$SD_markov_score[$q]{"g"}=log(($SD_markov_score[$q]{"g"}+0.000001)/$number_of_sample);	
	$SD_markov_score[$q]{"c"}=log(($SD_markov_score[$q]{"c"}+0.000001)/$number_of_sample);
    }

    for(my $key=0;$key<=($RBS_region_len-length($SD_seq));$key++){
	$loc_hash{$key}=log(($loc_hash{$key}+0.0001)/$number_of_sample);
    }
    
    #determine minimum markov score for threshold
    my $min_score=100000000;
    $max_score=-100000000000;
    for (my $n=0;$n<=$#RBS_seq;$n++){
	$tmp_score=&calculate_markov_score($RBS_seq[$n]);
	if ($tmp_score<$min_score){
	    $min_score=$tmp_score;
	}
	if ($tmp_score>$max_score){
	    $max_score=$tmp_score;
	}
    }

 


    my $min_loc_posib=100000000000;	
    my $max_loc_posib=-100000000000;	
    my $skip_val=log(0.0001/$number_of_sample);
    for(my $key=0;$key<=($RBS_region_len-length($SD_seq));$key++){
	next if ($loc_hash{$key} <= $skip_val); 
	if ($loc_hash{$key}>$max_loc_posib){
		$max_loc_posib=$loc_hash{$key};
	}
	if ($loc_hash{$key}<$min_loc_posib){
		$min_loc_posib=$loc_hash{$key};
	}
    }
    
     	
    $markov_score_threshold =($min_score+$max_score)*0.5+$min_loc_posib;
    $perfect_markov_score = $max_score+$max_loc_posib;
    $initial_threshold= ($perfect_markov_score+$markov_score_threshold)*0.5;

}
 
 
sub check_RBS_SD
{
    my $RBS_tmp=shift;
    my $threshold=shift;
    my $tmp_score=0;
    my $max_score=-100000;
    my $max_seq="";
    my $max_loc=0;	
    my $tmp_seq="";
    for (my $t=0;$t<=(length($RBS_tmp)-length($SD_seq));$t++){
	$tmp_seq=substr($RBS_tmp,$t,length($SD_seq));
	$tmp_score=&calculate_markov_score($tmp_seq)+$loc_hash{$t};
	if ($tmp_score>$max_score){
	    $max_seq=$tmp_seq;
	    $max_score=$tmp_score;
	    $max_loc=$t;	
	}
    }
    if ($max_score>$threshold){return ($max_seq." ".$max_loc." ".$max_score)}
    return "";
}

sub init_operons()
{
    open ORFLIST,$ARGV[1];
    my $p=0;
    my @orf_coords=();
    my @orf_name=();
    my @orf_reverse=();
    my $s;
    my $e;
 
    while (<ORFLIST>){
	chomp;
	next if ($_ eq "");
	next if ($_=~ /^\#/);
        #**** glimmer2 before get-putative
	if ($_ =~ /GC Proportion/){
	    while(1){
		$skip=<ORFLIST>;
		last if($skip=~ /Putative\sGenes/);
	    }
	    next;
	}
	#*********************************
	my @line_array=split(/\s+/,$_);
	if ($line_array[0] eq ""){
	        $s=$line_array[2];
		$e=$line_array[3];
		$orf_name[$p]=$line_array[1];
		$orf_reverse[$p]=($line_array[2]>$line_array[3])?1:0;	
	}
	else{
	    $s=$line_array[1];
	    $e=$line_array[2];
	    $orf_name[$p]=$line_array[0];
	    $orf_reverse[$p]=($line_array[1]>$line_array[2])?1:0;
	}	
	$orf_coords[$p][0]=($s>$e)?$e:$s;
	$orf_coords[$p++][1]=($s>$e)?$s:$e;
    }
 
#sort the coordinates according to the start coords of orfs
    my @sorted_orf_inds=sort {$orf_coords[$a][0] <=> $orf_coords[$b][0]} 0..$#orf_coords;
    
    my $start_operon=0;
#find if the orfs on forward strain are operons
    for (my $p=0;$p<=$#sorted_orf_inds;$p++){
	if ((!$start_operon) && (!$orf_reverse[$sorted_orf_inds[$p]])) {$start_operon=1; next;}
	if ($orf_reverse[$sorted_orf_inds[$p]]) {$start_operon=0; next;}
	if ($orf_reverse[$sorted_orf_inds[$p-1]]==$orf_reverse[$sorted_orf_inds[$p]]){
	    $is_operon{$orf_name[$sorted_orf_inds[$p]]}=10;
	}
	else{
	    $start_operon=0;
	}
}

#find if the orfs on reverse strain are operons
    my $start_operon=0;
    for (my $p=$#sorted_orf_inds;$p>=0;$p--){
	if ((!$start_operon) && ($orf_reverse[$sorted_orf_inds[$p]])) {$start_operon=1; next;}
	if (!$orf_reverse[$sorted_orf_inds[$p]]) {$start_operon=0; next;}
	if ($orf_reverse[$sorted_orf_inds[$p+1]]==$orf_reverse[$sorted_orf_inds[$p]]){
	    $is_operon{$orf_name[$sorted_orf_inds[$p]]}=10;
	}
	else{
	    $start_operon=0;
	}
    }

#initialize the boundry for the location to look for ,in finding start codon upstream
    for (my $p=0;$p<=$#sorted_orf_inds;$p++){
	if ($orf_reverse[$sorted_orf_inds[$p]]>0){ #if the current orf is reverse
	    if ($p == $#sorted_orf_inds){
		$orf_previous{$orf_name[$sorted_orf_inds[$p]]}=length($seq);
	    }
	    elsif ($orf_reverse[$sorted_orf_inds[$p+1]]>0){
		$orf_previous{$orf_name[$sorted_orf_inds[$p]]}=-1*($orf_coords[$sorted_orf_inds[($p+1)]][0]);
	    }
	    else{
		$orf_previous{$orf_name[$sorted_orf_inds[$p]]}=$orf_coords[$sorted_orf_inds[($p+1)]][0];
	    }
	}
	else{
	    if ($p == 0){
		$orf_previous{$orf_name[$sorted_orf_inds[$p]]}=1;
	    }
	    elsif ($orf_reverse[$sorted_orf_inds[$p-1]]>0){
		$orf_previous{$orf_name[$sorted_orf_inds[$p]]}=-1*$orf_coords[$sorted_orf_inds[$p-1]][1];
	    }
	    else{
		$orf_previous{$orf_name[$sorted_orf_inds[$p]]}=$orf_coords[$sorted_orf_inds[$p-1]][1];
	    }
	}   
    }
    close(ORFLIST);
}


sub calculate_threshold
{
    my $max_move=shift;
    my $moved=shift;
    return($markov_score_threshold+($perfect_markov_score-$markov_score_threshold)*($moved/$max_move));
}



sub sc_start{
	my $tmp_seq=shift;
	if ($tmp_seq eq "atg"){
		return 3;
	}
	if ($tmp_seq eq "gtg"){
		return 2;
	}
	if ($tmp_seq eq "ttg"){
		return 1;
	}
	if ($tmp_seq eq "cat"){
		return 3;
	}
	if ($tmp_seq eq "cac"){
		return 2;
	}
	if ($tmp_seq eq "caa"){
		return 1;
	}
	return 0;
}

