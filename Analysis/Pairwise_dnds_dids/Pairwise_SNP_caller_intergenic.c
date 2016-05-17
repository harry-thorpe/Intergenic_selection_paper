#include <stdio.h>
#include <string.h>
#include <math.h>

int main( int argc, char *argv[] )
{

// Collect command line args.
int arg_len = 1000;

char species[arg_len];
strcpy(species, argv[1]);

char analysis[arg_len];
strcpy(analysis, argv[2]);

char base_dir[arg_len];
strcpy(base_dir, argv[3]);

char log_file_name[arg_len];

// Open log file.
strcpy(log_file_name, base_dir);
strcat(log_file_name, "/Analysis/log.txt");

FILE *log_file = fopen (log_file_name, "a");

char output_file_name[arg_len];

// Open output file.
strcpy(output_file_name, base_dir);
strcat(output_file_name, "/Analysis/");
strcat(output_file_name, analysis);
strcat(output_file_name, "/");
strcat(output_file_name, species);
strcat(output_file_name, "_");
strcat(output_file_name, analysis);
strcat(output_file_name, "/");
strcat(output_file_name, species);
strcat(output_file_name, "_di.csv");

FILE *output = fopen (output_file_name, "w");

fprintf (output, "Isolate_1,Isolate_2,SNPs,Sites,Divergence\n");

// Open input isolate file.
char isolates_file_name[arg_len];

strcpy(isolates_file_name, base_dir);
strcat(isolates_file_name, "/Analysis/Core_genome_alignment/");
strcat(isolates_file_name, species);
strcat(isolates_file_name, "_Core_genome_alignment/");
strcat(isolates_file_name, species);
strcat(isolates_file_name, "_isolates.txt");

///media/harry/extra/Intergenic_variation_paper/Analysis/Core_genome_alignment/S_aureus_Core_genome_alignment/S_aureus_isolates.txt

// Count isolates.
char isolate_name[arg_len];

int isolate_count;
isolate_count = 0;
FILE *isolates_file = fopen (isolates_file_name, "r");
if(isolates_file != NULL){
	while(fgets(isolate_name, sizeof isolate_name, isolates_file) != NULL){
		isolate_count++;
	}
	fclose (isolates_file);
}else{
	perror(isolates_file_name);
}

// Make isolate name arrays.
char isolate_name_array[isolate_count][arg_len];
char isolate_file_name_array[isolate_count][arg_len];

// Populate isolate name arrays.
isolate_count = 0;
FILE *isolates_file_2 = fopen (isolates_file_name, "r");
if(isolates_file_2 != NULL){
	while(fgets(isolate_name, sizeof isolate_name, isolates_file_2) != NULL){
		int isolate_name_len;
		isolate_name_len = strlen(isolate_name);
		if(isolate_name[(isolate_name_len-1)] == '\n'){
			isolate_name[(isolate_name_len-1)] = '\0';
		}
		
		strcpy(isolate_name_array[isolate_count], isolate_name);
		
		strcpy(isolate_file_name_array[isolate_count], base_dir);
		strcat(isolate_file_name_array[isolate_count], "/Analysis/");
		strcat(isolate_file_name_array[isolate_count], analysis);
		strcat(isolate_file_name_array[isolate_count], "/");
		strcat(isolate_file_name_array[isolate_count], species);
		strcat(isolate_file_name_array[isolate_count], "_");
		strcat(isolate_file_name_array[isolate_count], analysis);
		strcat(isolate_file_name_array[isolate_count], "/Intergenic_files/");
		strcat(isolate_file_name_array[isolate_count], isolate_name);
		strcat(isolate_file_name_array[isolate_count], ".fasta");
		
		isolate_count++;
		
	}
	fclose (isolates_file_2);
}else{
	perror(isolates_file_name);
}

int a;
int b;
int i;
int genome_max = 1000000;
char genome1[genome_max];
char genome2[genome_max];

int count;
count = 0;
for(a=0; a<isolate_count; a++){
	FILE *file1 = fopen (isolate_file_name_array[a], "r");
	if(file1 != NULL){
		char line1[genome_max];
		while(fgets(line1, sizeof line1, file1) != NULL){
			if(line1[0] != '>'){
				int line1_len = strlen(line1);
				line1_len = (line1_len - 1);
				if(line1[line1_len] == '\n'){
					line1[line1_len] = '\0';
				}
			
				strcpy(genome1, line1);
			
				fclose (file1);
			}
		}
	}else{
		perror(isolate_file_name_array[a]);
	}
	for(b=0; b<isolate_count; b++){
		if(a > b){
			FILE *file2 = fopen (isolate_file_name_array[b], "r");
			if(file2 != NULL){
				char line2[genome_max];
				while(fgets(line2, sizeof line2, file2) != NULL){
					if(line2[0] != '>'){
						int line2_len = strlen(line2);
						line2_len = (line2_len - 1);
						if(line2[line2_len] == '\n'){
							line2[line2_len] = '\0';
						}
			
						strcpy(genome2, line2);
			
						fclose (file2);
					}
				}
			}else{
				perror(isolate_file_name_array[b]);
			}
			int genome_len;
			genome_len=strlen(genome1);
			
			int snp;
			snp = 0;
			
			int n;
			n = 0;
			
			for(i=0; i<genome_len; i++){
				if(genome1[i] == 'N' || genome2[i] == 'N'){
					n++;
				}else{
					if(genome1[i] != genome2[i]){
						snp++;
					}
				}
			}
			
			int sites;
			sites = (genome_len-n);
			
			float div;
			div = ((float) snp / sites);
			
			float JC_div;
			JC_div = (((float) -3 / 4) * (log(1 - (((float) 4 / 3) * div))));
			
			fprintf (output, "%s,%s,%i,%i,%0.6f\n", isolate_name_array[a], isolate_name_array[b], snp, sites, JC_div);
			
			count++;
			//if(count % 100 == 0){
			//	printf ("Pair %d completed.\n", count);
			//}
		}
	}
}

printf ("%s di calculated.\n", species);

fprintf (log_file, "%s di calculated.\n", species);

return 0;
}

