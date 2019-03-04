# ADARfootprints
multisequence alignments

* Takes multi nucleotide sequence fasta files as input and aligns properly framed and translated amino acid sequences and then back translates coding regions into nucleotide alignments.

* Takes advantage of publically available command line tools emboss and muscle tools.

___

## Getting Started
1.) Make sure you have emboss and muscle installed or you can use the provided VM



2.) Download tthe bash script.

3.) Make sure you have your .fasta nucleotide sequence file with a simple name in one folder and your corresponding .gff info file in a seperate folder.  You can run more then one multiple fasta file at time and it will keep them seperated. If you use the virtual machine these folders are on the desktop and once you have them ready run the following command.

```
bash /home/user/ADARfootprints.sh /home/user/Desktop/seq /home/user/Desktop/gff /home/user/ADARfootprints

```
   
   ***Argument 1 is the directory to find fasta files, argument 2 is the directory to find gff files, and argument 3 is the name of the directory you would like the results in.
   
___


## Output File

In ADARfootprints new directory you will find:

1.) Scripts with all the scripts and index files needed to run lit_search

2.) Categories folder csv file for graphs is the final folder

3.) Final folder with graphs for percent of total articles for each category compared between topics. A different bar graph for each topic.

4.) Topic_stats with the PMID for all the articles found for a topic.

5.) There will be a folder for each topic name from column 2 of your topics.csv file in the scripts folder.

   * In these folders you will find topicabPMID.csv which is the list of all the PMID that had abstracts with them,
   * topicabstracts.csv is the file with abstracts used for the search category search. 
   * ZIKVstat_names.csv is index file for data mining. 
   * ZIKVstats.csv which is metadata about each article downloaded.
   * There is also a folder called stats
      * In here is a folder for each meta data category country year accepted and year received (these can be changed manually in the Rscript)
      * In these folders you will find the pie chart and the csv files to make it.
   * There is also a folder called final
      * In here are the bargraphs for how many articles contained each search word by search category and topictotalwtot.tiff is the bargraph showing how many articles contained at least one search term shown by categories for the topic.
      * There is also a folder for each search category with a pie chart for how many articles contained at least one search term and the csv of all the unique PMID for that category.
      * Within charts folder you will find a pie chart for each individual search term number of articles compared to all for the topic.
      * Within final are the PMID list for each search term.
      * Within totals there are csv files for the charts in charts folder.

___

## Built With
EMBOSS - 


___

## References for tools.



