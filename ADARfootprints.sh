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
outdir="$3"
for gfffile in "$indirgff"/* ;
do
  f_name=$(echo "$gfffile" | cut -f 1 -d '.' | cut -f 6 -d'/')
  wkd="$outdir"/"$f_name"
  outdirtemp="$outdir"/temp
  for dirtomake in "$outdir" "$wkd" "$outdirtemp" ;
  do
    createdir
  done
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
    outdir="$3"
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
  for dirtomake in "$lwkd" "$twkd" "$fawkd" "$spwkd" "$CDSwkd" "$finald" ;
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
    move_file #renames
    file_in="$lwkd"/"$jname".fasta
    file_out="$lwkd"/"$jname"sixpack.txt
    tool=sixp
    run_tool #generate sixpack
    jname2=$(echo "$jname" | tr '[:upper:]' '[:lower:]')
    cat "$lwkd"/"$jname2".sixpack | awk 'f;/###/{f=1}' | sed '1,2d' | sed '$d' | sed '$d' | sed '$d' | sed 's/Total ORFs in frame //g' | sed 's/ :/,/g' | sed 's/   //g' | sed '$d' | sed 's/ //g' | sort -k2 -n -t, >> "$lwkd"/"$jname"frame.csv
    fnum=$(awk -F',' 'NR==1{print $1}' "$lwkd"/"$jname"frame.csv)
    file_in="$lwkd"/"$jname".fasta
    file_out="$twkd"/"$jname"trans.fasta
    tool=tran_seq
    run_tool #translates in the right frame according to sixpack
      addgaps #adds gaps to nucleotide sequences to shift to the right frame.
      rownum=$(awk -F',' '{if ($1 ~ /'$jname'/) print $4}' "$wkd"/"$f_name".csv)
      start=$(awk -F',' '{if ($1 ~ /'$jname'/) print $2}' "$wkd"/"$f_name".csv)
      end=$(awk -F',' '{if ($1 ~ /'$jname'/) print $3}' "$wkd"/"$f_name".csv)
      file_in="$lwkd"/"$jname".fasta
      file_out="$CDSwkd"/"$jname"CDS.fasta
      tool=tran_CDS #translates only CDS sequences
      run_tool
    mv "$lwkd"/"$jname"frame.csv "$spwkd"/
    mv "$lwkd"/"$jname2".sixpack "$spwkd"/
    mv "$lwkd"/"$jname".fasta "$fawkd"/
    mv "$lwkd"/"$jname"sixpack.txt "$spwkd"/
  done
  cat "$fawkd"/* >> "$wkd"/"$f_name"rawnuc.fasta #nuc with gaps for reading frame
  cat "$CDSwkd"/* >> "$wkd"/"$f_name"transCDS.fasta #translated CDS only
  cat "$twkd"/* >> "$wkd"/"$f_name"trans.fasta 
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
done
