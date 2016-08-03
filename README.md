# Intergenic selection paper

This repository contains the code for running the analysis in the paper.

##Dependencies

* Perl >= 5.8.3

##Instructions

###File structure

Clone the project into a local directory:

```bash
cd /somewhere/on/your/system
git clone https://github.com/harry-thorpe/Intergenic_selection_paper
```

Move into the project base directory:

```bash
cd Intergenic_selection_paper
```

Move into the Data directory and create input file directories:

```bash
cd Data
mkdir Reference_files
mkdir Alignment_files
mkdir GFF_files
mkdir Promoter_files
mkdir Terminator_files
cd ..
```

###Data files

To run the full analysis, you must transfer appropriate input files to the appropriate directories. The analysis can be run on multiple species, but the following instructions assume you are using `My_species`.

* Reference files - These should be in fasta format, with upper-case characters only (ATGCN), and the whole sequence on a single line. They should be called `My_species_reference.fasta`.
* Alignment files - These should be alignments in fasta format, aligned to the reference file, with upper-case characters only (ATGCN), and each sequence on a single line. The first sequence in the alignment should be the reference sequence. They should be called `My_species_alignment.fasta`.
* GFF files - These should be GFF files produced by Prokka. They should be called `My_species.gff`.
* Promoter files - These should be promoter annotation files for the reference genome used downloaded from http://pepper.molgenrug.nl/index.php/genome2d. They should be called `My_species_promoters.tab`.
* Terminator files - These should be terminator annotation files for the reference genome used downloaded from http://pepper.molgenrug.nl/index.php/genome2d. They should be called `My_species_terminators.tab`.

###Running the analysis

To run the analysis, you will have to edit two lines in the script `Intergenic_selection_paper/Analysis_script.sh`:

* Line 4 - Change the variable `$base_dir` to the directory where you have cloned the repository to.
* Line 7 - Change the species within `$species_array` to those which you want to analyse. These must be compatible with the input files.

Run the analysis:

```bash
bash Analysis_script.sh
```
