# ADARfootprints
multisequence alignments

* Takes multi nucleotide sequence fasta files as input and aligns properly framed and translated amino acid sequences and then back translates coding regions into nucleotide alignments.

* Takes advantage of publically available command line tools emboss and muscle tools.

___

## Getting Started
1.) Make sure you have emboss and muscle installed or you can use the provided VM



2.) Download tthe bash script.

3.) Make sure you have your .fasta nucleotide sequence file with a simple name in one folder and your corresponding .gff (with the same name) info file in a seperate folder.  You can run more then one multiple fasta file at time and it will keep them seperated. If you use the virtual machine these folders are on the desktop and once you have them ready run the following command.

```
bash /home/user/ADARfootprints.sh /home/user/Desktop/seq /home/user/Desktop/gff /home/user/ADARfootprints

```
   
   ***Argument 1 is the directory to find fasta files, argument 2 is the directory to find gff files, and argument 3 is the name of the directory you would like the results in.
   
___


## Output File

In ADARfootprints new directory you will find:

1.) final folder containing the nucleotide aligned sequences and the nucleotide CDS only aligned sequences for every set of fasta sequences provided in the seq input folder.

2.) there will be a folder for each set of fasta and gff files provided containing 

   * split folders containing individual sequence files for that set of fasta sequences provided.
   
      *In the splitfa folder, there are the raw sixplack files used to pick the best frame as well as translate folder containing the raw unaligned translated sequences used for the alignment 
   
   * a folder input with the raw input provided for that set of fasta sequences including a csv file generated by the tool for determining CDS.
   
   * a folder raw_seq which contains intermediate files
   
      *files with trans in the name indicate translated sequences
      *nuc indicate nucleotide sequences 
      *backtran indicate protein alignment back tranlated to nucleotide
      *CDS indicate that the sequences have been filtered to only contain CDS sequences according to the csv file generated.

___

## Built With

EMBOSS: http://emboss.open-bio.org/
Muscle: http://www.drive5.com/muscle/muscle.html


___

## References for tools.

Rice P., Bleasby A and Ison J. The EMBOSS Users Guide. Cambridge University Press
Ison J., Rice P. and Bleasby A. The EMBOSS Developers Guide. Cambridge University Press

Bleasby A., Ison J. and Rice P. The EMBOSS Administrators Guide. Cambridge University Press

Rice P., Longden I. and Bleasby A. EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics. 2000 16(6):276-277

Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97.
