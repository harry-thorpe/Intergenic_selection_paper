#!/usr/bin/bash

species=$1
analysis=$2
base_dir=$3

half=$4

species_analysis="$species""_""$analysis"

IFS=$'\n' read -d '' -r -a isolate_array < "${base_dir}/Analysis/Core_genome_alignment/${species}""_Core_genome_alignment/${species}""_isolates.txt"

isolate_count=${#isolate_array[@]}

echo $isolate_count

mkdir "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half"

cp "${base_dir}/Analysis/${analysis}/yn00.ctl" "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half"

cd "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half"

for ((a=0; a<$isolate_count; a++));
do
	for ((b=0; b<$isolate_count; b++));
	do
		if [ $a -gt $b ];
		then
			iso_1=${isolate_array[$a]}
			iso_2=${isolate_array[$b]}
			
			cat "${base_dir}/Analysis/${analysis}/${species_analysis}/Gene_files_half_$half/$iso_1.fasta" >> "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half/ali.fasta"
			cat "${base_dir}/Analysis/${analysis}/${species_analysis}/Gene_files_half_$half/$iso_2.fasta" >> "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half/ali.fasta"
	
			yn00 yn00.ctl
			
			iso_2_regex="^""$iso_2"" \+[^ ]\+"
			
			res=$(grep "$iso_2_regex" "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half/dnds_out.txt")
	
			echo "$iso_1            $res" >> "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half/${species}""_dnds_half_$half.txt"
	
			rm "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half/ali.fasta"
			
		fi
	done
done

cp "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half/${species}""_dnds_half_$half.txt" "${base_dir}/Analysis/${analysis}/${species_analysis}/${species}""_dnds_half_$half.txt"

rm -r "${base_dir}/Analysis/${analysis}/${species_analysis}/dnds_tmp_$half"

cd "${base_dir}/Analysis/${analysis}"

echo "$species dnds calculated."

echo "$species dnds calculated." >> "$base_dir/Analysis/log.txt"

