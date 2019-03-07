#!/usr/bin/env bash

##need to enter directory where input fasta are located first, then where gff are and then second argument is the directory of where to put the output files.
createdir() {
if [ ! -d "$dirtomake" ];
then
  mkdir "$dirtomake"
fi
}
move_file() {
if [ -s "$file_in" ];
then
  if [ ! -s "$file_out" ];
  then
    mv "$file_in" "$file_out"
  fi
fi
}
remove_file() {
if [ -s "$file_out" ];
then
  rm "$file_in"
fi
}
mus_align() {
muscle -group -quiet -in "$file_in" -out "$file_out"
}
tran_seq() {
if [ "$fnum" == "4" ];
then
  fnum=$(echo "-1")
fi
if [ "$fnum" == "5" ];
then
  fnum=$(echo "-2")
fi
if [ "$fnum" == "6" ];
then
  fnum=$(echo "-3")
fi
transeq -auto -sequence "$file_in" -outseq "$file_out" -frame="$fnum"
}
tran_CDS() {
transeq -auto -sequence "$file_in" -outseq "$file_out" -reg="$start"-"$end"
}
back_tran() {
backtranseq -auto -sequence "$file_in" -outfile "$file_out"
}
sixp() {
sixpack -auto -sequence "$file_in" -outseq "$file_out"
}
tran_align() {
tranalign -auto -asequence "$file_in" -bsequence "$file_in2" -outseq "$file_out" 
}
run_tool() {
if [ ! -s "$file_out" ];
then
  if [ -s "$file_in" ];
  then
    "$tool"
  fi
fi
}
run_tool2() {
if [ ! -s "$file_out" ];
then
  if [[ -s "$file_in" && -s "$file_in2" ]];
  then
    "$tool"
  fi
fi
}
addgaps() {
if [ "$fnum" == "3" ];
    then
      cat "$lwkd"/"$jname".fasta | sed '2s/.*/-&/' >> "$lwkd"/temp.fasta
      rm "$lwkd"/"$jname".fasta
      mv "$lwkd"/temp.fasta "$lwkd"/"$jname".fasta
    fi
    if [ "$fnum" == "2" ];
    then
      cat "$lwkd"/"$jname".fasta | sed '2s/.*/--&/' >> "$lwkd"/temp.fasta
      rm "$lwkd"/"$jname".fasta
      mv "$lwkd"/temp.fasta "$lwkd"/"$jname".fasta
    fi
    if [[ "$fnum" == "4" || "$fnum" == "5" || "$fnum" == "6" ]];
    then
      revseq "$lwkd"/"$jname".fasta -sequence "$lwkd"/"$jname".fasta -outseq "$lwkd"/"$jname"rev.fasta -nocomp 
      rm "$lwkd"/"$jname".fasta
      mv "$lwkd"/"$jname"rev.fasta "$lwkd"/"$jname".fasta    
      if [ "$fnum" == "6" ];
      then
        cat "$lwkd"/"$jname".fasta | sed '2s/.*/-&/' >> "$lwkd"/temp.fasta
        rm "$lwkd"/"$jname".fasta
        mv "$lwkd"/temp.fasta "$lwkd"/"$jname".fasta
      fi
      if [ "$fnum" == "5" ];
      then
        cat "$lwkd"/"$jname".fasta | sed '2s/.*/--&/' >> "$lwkd"/temp.fasta
        rm "$lwkd"/"$jname".fasta
        mv "$lwkd"/temp.fasta "$lwkd"/"$jname".fasta
      fi      
    fi
}
indirseq="$1"
indirgff="$2"
indirmeta="$3"
outdir="$4"
#For country stats column_num=10
#For vector column_num=8
#For host column_num=9
#For collection data=7
#For sequence length=6
#For genbanksubtype=4
#For VIPRsubtype=3
################################################################################################################################################
#this summaries metadata about sequences can be used later to add filtering out a subseq of sequences to work with
################################################################################################################################################
for metafile in "$indirmeta"/* ; 
do
  f_name=$(echo "$metafile" | cut -f 1 -d '.' | cut -f 6 -d'/')
  metawkd="$outdir"/metadata
  for dirtomake in "$outdir" "$metawkd" ;
  do
    createdir
  done
  file_in="$indirmeta"/"$f_name".tsv
  file_out="$outdir"/metadata/"$f_name".tsv
  cp "$file_in" "$file_out"
  cat "$outdir"/metadata/"$f_name".tsv | sed 's/ /_/g' | cut -d$'\t' -f5,3,4,6,7,8,9,10 | sed 's/\t/,/g' | sed 's/*/_/g' | sed 's/-//g' | sed 's/\//_/g' >> "$outdir"/metadata/"$f_name"meta.csv
  for col_num in 2 3 4 5 6 7 8 ;
  do
    echo ""$col_num""
    metaname=$(awk 'NR==1{print $'$col_num'}' "$outdir"/metadata/"$f_name"meta.csv) 
    echo ""$metaname""
    cat "$outdir"/metadata/"$f_name"meta.csv | cut -d',' -f$col_num | sort | uniq -ci | sed 's/ \+/,/g' | sed '1i name,freq' > "$outdir"/metadata/"$f_name""$col_num".csv
    #file_in=
    #file_out=
    #Rscript barchart.R "$file_in" "$file_out"
  done
done
################################################################################################################################################
#this makes the file to filter sequences to just CDS( only one for now) maybe later add option to have multiple ORFs for not flavi
################################################################################################################################################
for gfffile in "$indirgff"/* ; 
do
  tempdir="$outdir"/temp
  wkd="$outdir"/"$f_name"
  for dirtomake in "$tempdir" "$wkd" ;
  do
    createdir
  done
  f_name=$(echo "$gfffile" | cut -f 1 -d '.' | cut -f 6 -d'/')
  wkd="$outdir"/"$f_name"
  outdirtemp="$outdir"/temp
  file_in="$indirgff"/"$f_name".gff
  file_out="$wkd"/"$f_name".gff
  cp "$indirgff"/"$f_name".gff "$wkd"/"$f_name".gff
  cat "$wkd"/"$f_name".gff | sed '1,3d' | sed '/^##s/ d' | sed '/##FASTA/,$d' >> "$outdir"/temp/"$f_name"temp.txt
  cat "$outdir"/temp/"$f_name"temp.txt | awk '/ViPR	CDS/{print NR}' >> "$outdir"/temp/"$f_name"temp2.csv
  cat "$outdir"/temp/"$f_name"temp2.csv | awk '{$2=$1-1}1' | sed 's/ /,/g' >> "$outdir"/temp/"$f_name"temp3.csv
  INPUT="$outdir"/temp/"$f_name"temp3.csv  
  {
  [ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
  while IFS=, read -r CDS ID
  do
    outdir="$4"
    line1=$(sed "${ID}q;d" "$outdir"/temp/"$f_name"temp.txt)
    line2=$(sed "${CDS}q;d" "$outdir"/temp/"$f_name"temp.txt)
    echo ""$line1"" >> "$outdir"/temp/"$f_name"tempID.csv
    echo ""$line2"" >> "$outdir"/temp/"$f_name"tempCDS.csv
  done
  } < $INPUT
  IFS=$OLDIFS
  mv "$outdir"/temp/"$f_name"tempID.csv "$wkd"/"$f_name"tempID.csv
  mv "$outdir"/temp/"$f_name"tempCDS.csv "$wkd"/"$f_name"tempCDS.csv
  cat "$wkd"/"$f_name"tempID.csv | sed 's/ /,/g' | sed 's/	/,/g' | sed 's/  /,/g' | cut -d ',' -f1 | sed 's/$/,/g' >> "$wkd"/"$f_name"tempIDonly.csv
  cat "$wkd"/"$f_name"tempCDS.csv | sed 's/ /,/g' | sed 's/	/,/g' | sed 's/  /,/g' | cut -d ',' -f4,5 >> "$wkd"/"$f_name"tempCDSonly.csv
  paste -d, "$wkd"/"$f_name"tempIDonly.csv <(cut -d, -f1,2- "$wkd"/"$f_name"tempCDSonly.csv) | awk -F ',' '{$5=NR}1' | cut -d ',' -f5,1,3,4 | sed 's/  /,/g' | sed 's/ /,/g' >> "$wkd"/"$f_name".csv
done
################################################################################################################################################
#this runs prep for alignments: frame shift, translation, and CDS filtering
################################################################################################################################################
for seqfile in "$indirseq"/* ;
do
  f_name=$(echo "$seqfile" | cut -f 1 -d '.' | cut -f 6 -d'/')
  echo ""$f_name""
  wkd="$outdir"/"$f_name"
  lwkd="$wkd"/splitfa
  twkd="$lwkd"/translate
  fawkd="$lwkd"/fasta_files
  spwkd="$lwkd"/sixpack_files
  CDSwkd="$wkd"/splitCDS
  finald="$outdir"/final
  inputd="$wkd"/input
  rawd="$wkd"/raw_seq
  for dirtomake in "$lwkd" "$twkd" "$fawkd" "$spwkd" "$CDSwkd" "$finald" "$inputd" "$rawd" ;
  do
    createdir
  done
  file_in="$indirseq"/"$f_name".fasta
  file_out="$wkd"/"$f_name".fasta
  cp "$file_in" "$file_out" 
  cd "$lwkd"
  csplit -s -z "$wkd"/"$f_name".fasta /\>.*/ {*} #splits file
  for fa_file in xx* ; 
  do
    name=$(head -1 $fa_file)
    jname=$(echo $name | cut -d'>' -f 2)
    file_in="$fa_file"
    file_out="$lwkd"/"$jname".fasta
    move_file #renames file to accension number
    file_in="$lwkd"/"$jname".fasta
    file_out="$lwkd"/"$jname"sixpack.txt
    tool=sixp
    run_tool #generate sixpack for each sequence
    jname2=$(echo "$jname" | tr '[:upper:]' '[:lower:]')
    cat "$lwkd"/"$jname2".sixpack | awk 'f;/###/{f=1}' | sed '1,2d' | sed '$d' | sed '$d' | sed '$d' | sed 's/Total ORFs in frame //g' | sed 's/ :/,/g' | sed 's/   //g' | sed '$d' | sed 's/ //g' | sort -k2 -n -t, >> "$lwkd"/"$jname"frame.csv #takes sixpack picks which frame is the best based on least number of stop codons
    fnum=$(awk -F',' 'NR==1{print $1}' "$lwkd"/"$jname"frame.csv)
    file_in="$lwkd"/"$jname".fasta
    file_out="$twkd"/"$jname"trans.fasta
    tool=tran_seq
    run_tool #translates in the right frame according to lowest number of stop codons in the sixpack
      addgaps #adds gaps to nucleotide sequences to shift to the right frame for pulling CDS sequences later
      rownum=$(awk -F',' '{if ($1 ~ /'$jname'/) print $4}' "$wkd"/"$f_name".csv)
      start=$(awk -F',' '{if ($1 ~ /'$jname'/) print $2}' "$wkd"/"$f_name".csv)
      end=$(awk -F',' '{if ($1 ~ /'$jname'/) print $3}' "$wkd"/"$f_name".csv)
      file_in="$lwkd"/"$jname".fasta
      file_out="$CDSwkd"/"$jname"CDS.fasta
      tool=tran_CDS 
      run_tool #translates only CDS sequences to align later
    mv "$lwkd"/"$jname"frame.csv "$spwkd"/
    mv "$lwkd"/"$jname2".sixpack "$spwkd"/
    mv "$lwkd"/"$jname".fasta "$fawkd"/
    mv "$lwkd"/"$jname"sixpack.txt "$spwkd"/
  done
################################################################################################################################################
#aligns protein sequences and back translates to nucleotide alignments
################################################################################################################################################
  cat "$lwkd"/* >> "$wkd"/"$f_name"rawnuc.fasta #create big nuc with gaps for reading frame file
  cat "$CDSwkd"/* >> "$wkd"/"$f_name"transCDS.fasta #create big translated CDS only file
  cat "$twkd"/* >> "$wkd"/"$f_name"trans.fasta #create big translated whole sequence file
  file_in="$wkd"/"$f_name"trans.fasta
  file_out="$wkd"/"$f_name"transalign.fasta
  tool=mus_align # aligns translated sequences
  run_tool
  file_in="$wkd"/"$f_name"transalign.fasta
  file_out="$wkd"/"$f_name"transalignbacktran.fasta
  tool=back_tran # back translates to nucleotides
  run_tool
  file_in="$wkd"/"$f_name"transalignbacktran.fasta
  file_in2="$wkd"/"$f_name"transalign.fasta
  file_out="$outdir"/final/"$f_name"nucalign.fasta
  tool=tran_align # aligns nucleotides based on amino acid alignment
  run_tool2
  file_in="$wkd"/"$f_name"transCDS.fasta
  file_out="$wkd"/"$f_name"transCDSalign.fasta
  tool=mus_align # alignes CDS only translated sequences
  run_tool
  file_in="$wkd"/"$f_name"transCDSalign.fasta
  file_out="$wkd"/"$f_name"transCDSalignbacktran.fasta
  tool=back_tran # back translates to nucleotides 
  run_tool
  file_in="$wkd"/"$f_name"transCDSalignbacktran.fasta
  file_in2="$wkd"/"$f_name"transCDSalign.fasta
  file_out="$outdir"/final/"$f_name"nucalignCDS.fasta
  tool=tran_align # aligns nucleotides based on amino acid alignment for CDS only
  run_tool2
  cd
################################################################################################################################################
#clean up
################################################################################################################################################
  rm "$wkd"/"$f_name"temp*
  mv "$wkd"/"$f_name".*  "$wkd"/input/
  mv "$wkd"/"$f_name"* "$wkd"/raw_seq/
  cat "$outdir"/final/*filteredfasta.csv | sed '1s Alignment,ID,error_type' >> "$outdir"/final/filteredfasta.csv
  virus=$(echo ""$f_name"")
  sixpack=$(grep -c "sixpack" "$outdir"/final/"$f_name"filteredfasta.csv)
  CDS=$(grep -c "CDS" "$outdir"/final/"$f_name"filteredfasta.csv)
  unknown=$(grep -c "unknown" "$outdir"/final/"$f_name"filteredfasta.csv)
  totalstart=$(grep -c ">" "$wkd"/input/"$f_name".fasta)
  totalfiltered=$(wc -l "$outdir"/final/"$f_name"filteredfasta.csv)
  #totalused=$(ls -1 | wc -l "$fawkd"/*) file_count=$( shopt -s nullglob ; set -- $directory_to_search_inside/* ; echo $#)
  #percentseq=
  #numofsites=
  #variablesites=
  if [ ! -f "$outdir"/final/alignmentsummary.csv ];
  then
    echo "virus_name,sixpack_error,CDSnotfound,toomanyunknownAA,seqdownloaded,seqfiltered,seqused,percentused,totalsitesalignment,variable_sites" >> "$outdir"/final/alignmentsummary.csv
  fi
  echo ""$virus","$sixpack","$CDS","$unknown","$totalstart","$totalfiltered","$totalused","$percentseq","$numofsites","$variablesites"" >> "$outdir"/final/alignmentsummary.csv

done
rm -rf "$outdir"/temp
