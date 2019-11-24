#!/usr/bin/env bash

##need to enter directory where input fasta are located first, then where gff are and then second argument is the directory of where to put the output files.
checktools() {
if ! [ -x "$(command -v muscle)" ]; then
  echo 'Error: muscle is not installed.' >&2
  sudo apt-get install muscle
fi
if ! [ -x "$(command -v transeq)" ]; then
  echo 'Error: emboss is not installed.' >&2
  sudo apt-get install emboss
fi
if ! [ -x "$(command -v csvtool)" ]; then
  echo 'Error: emboss is not installed.' >&2
  sudo apt-get install csvtool
fi
}
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
echo "You need to have files already downloads from www.viprbrc.org
You need to download sequences with only GenBank Accession no other information in the fasta file format.
You need to download gff file indicating the CDS.
You need to download metadata in the tsv file format.
Are all of these files in place? Please choose one of the following:
yes or no"
read answer
if [ "$answer" == "yes" ];
then
  echo "Where are your input files located?
fasta format sequences only annotated with GenBank Accession should be in /path/to/files/sequences
gff files should be located in /path/to/files/gff
metadata should be in tsv format in /path/to/files/metadata
all three of these folders need to be in the directory you enter here. For example /home/user/raw_data"
  read indir
  echo "Where would you like the results created? Please enter a directory /path/to/Results"
  read outdir
else
  echo "Please download files and try to run ADARfootprints when all the files are ready."
  exit
fi
checktools
indirseq="$indir"/sequences
indirgff="$indir"/gff
indirmeta="$indir"/metadata
if [ -s ~/ADARfootprints/config.cfg ];
then
  rm ~/ADARfootprints/config.cfg
fi
if [ -s ~/ADARfootprints/config.cfg.defaults ];
then
  rm ~/ADARfootprints/config.cfg.defaults
fi
echo "indir="$indir"
outdir="$outdir"
indirseq="$indir"/sequences
indirgff="$indir"/gff
indirmeta="$indir"/metadata" >> ~/ADARfootprints/config.cfg
echo "indirseq=Default Value
indirgff=Default Value
indirmeta=Default Value" >> ~/ADARfootprints/config.cfg.defaults
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
  file_in="$metafile"
  namefile=$(echo "${file_in##*/}")
  f_name=$(echo "${namefile%%.*}")
  #extension=$(echo "${bamfiles##*.}")
  #f_name=$(echo "$metafile" | cut -f 1 -d '.' | cut -f 6 -d'/')
  metawkd="$outdir"/metadata
  dirtomake="$outdir"
  createdir
  dirtomake="$metawkd"
  createdir
  file_in="$indirmeta"/"$f_name".tsv
  file_out="$outdir"/metadata/"$f_name".tsv
  cp "$file_in" "$file_out"
  cat "$outdir"/metadata/"$f_name".tsv | sed 's/ /_/g' | cut -d$'\t' -f5,3,4,6,7,8,9,10 | sed 's/\t/,/g' | sed 's/*/_/g' | sed 's/-//g' | sed 's/\//_/g' >> "$outdir"/metadata/"$f_name"meta.csv
  for col_num in 2 3 4 5 6 7 8 ;
  do
    echo ""$col_num""
    metaname=$(awk 'NR==1{print $'$col_num'}' "$outdir"/metadata/"$f_name"meta.csv) 
    echo ""$metaname""
    cat "$outdir"/metadata/"$f_name"meta.csv | cut -d',' -f$col_num | uniq -ci | sed 's/ \+/,/g' | sed '1i name,freq' > "$outdir"/metadata/"$f_name""$col_num".csv
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
  file_in="$gfffile"
  namefile=$(echo "${file_in##*/}")
  f_name=$(echo "${namefile%%.*}")
  tempdir="$outdir"/temp
  dirtomake="$tempdir"
  createdir
  wkd="$outdir"/"$f_name"
  dirtomake="$wkd"
  createdir
  #echo "f_name="$f_name"" >> ~/ADARfootprints/config.cfg.defaults
  #extension=$(echo "${bamfiles##*.}")
  #f_name=$(echo "$gfffile" | cut -f 1 -d '.' | cut -f 6 -d'/')
  file_in="$indirgff"/"$f_name".gff
  file_out="$wkd"/"$f_name".gff
  cp "$indirgff"/"$f_name".gff "$wkd"/"$f_name".gff
  cat "$wkd"/"$f_name".gff | sed '1,3d' | sed '/^##s/ d' | sed '/##FASTA/,$d' >> "$outdir"/temp/"$f_name"temp.txt
  cat "$outdir"/temp/"$f_name"temp.txt | awk '/ViPR	CDS/{print NR}' >> "$outdir"/temp/"$f_name"temp2.csv
  cat "$outdir"/temp/"$f_name"temp2.csv | awk '{$2=$1-1}1' | sed 's/ /,/g' >> "$outdir"/temp/"$f_name"temp3.csv
  #cd ~/ADARfootprints
  INPUT="$outdir"/temp/"$f_name"temp3.csv  
  {
  [ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
  while IFS=, read -r CDS ID
  do
   # source config.shlib; # load the config library functions
   # outdir="$(config_get outdir)";
   # f_name="$(config_get f_name)";
    outdir=$outdir
    f_name=$f_name
    line1=$(sed "${ID}q;d" "$outdir"/temp/"$f_name"temp.txt)
    line2=$(sed "${CDS}q;d" "$outdir"/temp/"$f_name"temp.txt)
    echo ""$line1"" >> "$outdir"/temp/"$f_name"tempID.csv
    echo ""$line2"" >> "$outdir"/temp/"$f_name"tempCDS.csv
  done
  } < $INPUT
  IFS=$OLDIFS
  #cat ~/ADARfootprints/config.cfg.defaults | sed '/f_name=/d' >> ~/ADARfootprints/temp.cfg
  #if [ -s ~/ADARfootprints/temp.cfg ];
  #then
  #  rm ~/ADARfootprints/config.cfg.defaults
  #  mv ~/ADARfootprints/temp.cfg ~/ADARfootprints/config.cfg.defaults
  #fi
  mv "$outdir"/temp/"$f_name"tempID.csv "$wkd"/"$f_name"tempID.csv
  mv "$outdir"/temp/"$f_name"tempCDS.csv "$wkd"/"$f_name"tempCDS.csv
  cat "$wkd"/"$f_name"tempID.csv | sed 's/ /,/g' | sed 's/	/,/g' | sed 's/  /,/g' | cut -d ',' -f1 | sed 's/$/,/g' >> "$wkd"/"$f_name"tempIDonly.csv
  cat "$wkd"/"$f_name"tempCDS.csv | sed 's/ /,/g' | sed 's/	/,/g' | sed 's/  /,/g' | cut -d ',' -f4,5 >> "$wkd"/"$f_name"tempCDSonly.csv
  paste -d, "$wkd"/"$f_name"tempIDonly.csv <(cut -d, -f1,2- "$wkd"/"$f_name"tempCDSonly.csv) | awk -F ',' '{$5=NR}1' | cut -d ',' -f5,1,3,4 | sed 's/  /,/g' | sed 's/ /,/g' | awk -F',' 'BEGIN { OFS = "," } { print $1,$2,$3,$4,$3-$2 } 1' | awk -F',' 'NR%2==1' >> "$wkd"/"$f_name".csv
  if [ "$ORF" == "2" ];
  then
    cat "$wkd"/"$f_name".csv | awk -F',' 'NR%2==1' >> "$f_name"ORF1.csv
    cat "$wkd"/"$f_name".csv | awk -F',' 'NR%2==2' >> "$f_name"ORF2.csv
  fi
done
################################################################################################################################################
#this runs prep for alignments: frame shift, translation, and CDS filtering
################################################################################################################################################
for seqfile in "$indirseq"/* ;
do
  file_in="$seqfile"
  namefile=$(echo "${file_in##*/}")
  f_name=$(echo "${namefile%%.*}")
  #f_name=$(echo "$seqfile" | cut -f 1 -d '.' | cut -f 6 -d'/')
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
    dirtomake="$dirtomake"
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
    if [ "$ORF" == "2" ];
    then
      rownum=$(awk -F',' '{if ($1 ~ /'$jname'/) print $4}' "$wkd"/"$f_name"ORF1.csv)
      start=$(awk -F',' '{if ($1 ~ /'$jname'/) print $2}' "$wkd"/"$f_name"ORF1.csv)
      end=$(awk -F',' '{if ($1 ~ /'$jname'/) print $3}' "$wkd"/"$f_name"ORF1.csv)
      file_in="$lwkd"/"$jname".fasta
      file_out="$CDSwkd"/"$jname"CDSORF1.fasta
      tool=tran_CDS 
      run_tool #translates only CDS sequences to align later
      rownum=$(awk -F',' '{if ($1 ~ /'$jname'/) print $4}' "$wkd"/"$f_name"ORF2.csv)
      start=$(awk -F',' '{if ($1 ~ /'$jname'/) print $2}' "$wkd"/"$f_name"ORF2.csv)
      end=$(awk -F',' '{if ($1 ~ /'$jname'/) print $3}' "$wkd"/"$f_name"ORF2.csv)
      file_in="$lwkd"/"$jname".fasta
      file_out="$CDSwkd"/"$jname"CDSORF2.fasta
      tool=tran_CDS 
      run_tool #translates only CDS sequences to align later
    else
      rownum=$(awk -F',' '{if ($1 ~ /'$jname'/) print $4}' "$wkd"/"$f_name".csv)
      start=$(awk -F',' '{if ($1 ~ /'$jname'/) print $2}' "$wkd"/"$f_name".csv)
      end=$(awk -F',' '{if ($1 ~ /'$jname'/) print $3}' "$wkd"/"$f_name".csv)
      file_in="$lwkd"/"$jname".fasta
      file_out="$CDSwkd"/"$jname"CDS.fasta
      tool=tran_CDS 
      run_tool #translates only CDS sequences to align later
    fi
    mv "$lwkd"/"$jname"frame.csv "$spwkd"/
    mv "$lwkd"/"$jname2".sixpack "$spwkd"/
    mv "$lwkd"/"$jname".fasta "$fawkd"/
    mv "$lwkd"/"$jname"sixpack.txt "$spwkd"/
  done
################################################################################################################################################
#aligns protein sequences and back translates to nucleotide alignments
################################################################################################################################################
  cat "$lwkd"/* >> "$wkd"/"$f_name"rawnuc.fasta #create big nuc with gaps for reading frame file
  if [ "$ORF" == "2" ];
  then
    cat "$CDSwkd"/* >> "$wkd"/"$f_name"transCDSORF1.fasta #create big translated CDS only file
    file_in="$wkd"/"$f_name"transCDSORF1.fasta
    file_out="$wkd"/"$f_name"transCDSalignORF1.fasta
    tool=mus_align # alignes CDS only translated sequences
    run_tool
    file_in="$wkd"/"$f_name"transCDSalignORF1.fasta
    file_out="$wkd"/"$f_name"transCDSalignbacktranORF1.fasta
    tool=back_tran # back translates to nucleotides 
    run_tool
    file_in="$wkd"/"$f_name"transCDSalignbacktranORF1.fasta
    file_in2="$wkd"/"$f_name"transCDSalignORF1.fasta
    file_out="$outdir"/final/"$f_name"nucalignCDSORF1.fasta
    tool=tran_align # aligns nucleotides based on amino acid alignment for CDS only
    run_tool2
    cat "$CDSwkd"/* >> "$wkd"/"$f_name"transCDSORF2.fasta #create big translated CDS only file
    file_in="$wkd"/"$f_name"transCDSORF2.fasta
    file_out="$wkd"/"$f_name"transCDSalignORF2.fasta
    tool=mus_align # alignes CDS only translated sequences
    run_tool
    file_in="$wkd"/"$f_name"transCDSalignORF2.fasta
    file_out="$wkd"/"$f_name"transCDSalignbacktranORF2.fasta
    tool=back_tran # back translates to nucleotides 
    run_tool
    file_in="$wkd"/"$f_name"transCDSalignbacktranORF2.fasta
    file_in2="$wkd"/"$f_name"transCDSalignORF2.fasta
    file_out="$outdir"/final/"$f_name"nucalignCDSORF2.fasta
    tool=tran_align # aligns nucleotides based on amino acid alignment for CDS only
    run_tool2
  else 
    for files in "$CDSwkd"/* ;
    do
      file_in="$files"
      namefile=$(echo "${file_in##*/}")
      seq_name=$(echo "${namefile%%.*}")
      seq_nameshort=${seq_name::-3}
      missing=$(awk -F',' '$1 == "'$seq_nameshort'" {print $0}' "$wkd"/"$f_name".csv)
      if [ "$missing" == "" ];
      then
        rm "$CDSwkd"/"$seq_name".fasta
        echo "$seq_nameshort,removed" >> "$wkd"/"$f_name"removedCDS.csv
      else
        echo "$seq_nameshort,kept" >> "$wkd"/"$f_name"removedCDS.csv
      fi
      missingAA=$(awk -F',' '/XXX/{print$0}' "$CDSwkd"/"$seq_name".fasta)
      echo "$missingAA"
      if [ "$missing" == "" ];
      then
        rm "$CDSwkd"/"$seq_name".fasta
        echo "$seq_nameshort,removed" >> "$wkd"/"$f_name"removedmissingAA.csv
      else
        echo "$seq_nameshort,kept" >> "$wkd"/"$f_name"removedmissingAA.csv        
      fi
      toomanystops=$(awk -F',' 'NR==1{print$2}' "$spwkd"/"$seq_nameshort"frame.csv)
      toomanystops=$(expr "$toomanystops" - 1)
      if [ "$toomanystops" -ge "20" ];
      then
        rm "$CDSwkd"/"$seq_name".fasta
        echo "$seq_nameshort,removed" >> "$wkd"/"$f_name"removedtoomanystops.csv
      else
        echo "$seq_nameshort,kept" >> "$wkd"/"$f_name"removedtoomanystops.csv        
      fi
    done
    cat "$CDSwkd"/* >> "$wkd"/"$f_name"transCDS.fasta #create big translated CDS only file
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
  fi
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
  cd
################################################################################################################################################
#transform fasta to matrix to add ADAR status then back to fasta for DAMBE input
################################################################################################################################################
  file_in="$outdir"/final/"$f_name"nucalignCDS.fasta
  file_out="$outdir"/final/"$f_name"nucalignCDSlinear.fasta
  while read line;do if [ "${line:0:1}" == ">" ]; then echo -e "\n"$line; else echo $line | tr -d '\n' ; fi; done < "$file_in" > "$file_out"
  file_in="$outdir"/final/"$f_name"nucalignCDSlinear.fasta
  file_out="$outdir"/final/"$f_name"nucalignCDS.csv
  cat "$file_in" | sed '1d' | sed 's/$/,/' | sed 'N;s/\n/ /' >> "$file_out"
  file_in="$outdir"/final/"$f_name"nucalignCDS.csv
  file_out="$outdir"/final/"$f_name"nucalignCDSRownames.csv
  file_out2="$outdir"/final/"$f_name"nucalignCDSNoRownames.csv
  file_out3="$outdir"/final/"$f_name"nucalignCDSfinal.csv
  cat "$file_in" | sed 's/,[^,]*,/,/' >> "$file_out"
  cat "$file_in" | sed 's/[^,]*,//' | sed 's/A/A,/g' | sed 's/C/C,/g' | sed 's/G/G,/g' | sed 's/T/T,/g' | sed 's/-/-,/g' | sed 's/N/N,/g' >> "$file_out2"
paste -d' ' "$file_out" "$file_out2" >> "$file_out3"
  if [ -s "$file_out3" ];
  then
    rm "$file_out"
    rm "$file_out2"
  fi
  file_in="$outdir"/final/"$f_name"nucalignCDSfinal.csv
  file_out="$outdir"/final/"$f_name"nucalignCDStrans.csv
  file_out2="$outdir"/final/"$f_name"nucalignCDStrans2.csv
  cat "$file_in" | csvtool transpose - >> "$file_out"
  header=$(head -n 1 "$file_out")
  cat "$file_out" | awk 'NR > 1{printf "%s,%s\n", NR,$0}' | awk -F "," '{if($2=="A")print $0",1";else print $0",0"}' | sed '1i site,'$header',ADAReditingsite' >> "$file_out2"



################################################################################################################################################
#clean up
################################################################################################################################################
  rm "$wkd"/"$f_name"temp*
  mv "$wkd"/"$f_name".*  "$wkd"/input/
  mv "$wkd"/"$f_name"* "$wkd"/raw_seq/
  file_in="$outdir"/final/"$f_name"nucalignCDSfinal.csv
  file_out="$outdir"/final/"$f_name"nucalignCDStrans.csv
  file_out2="$outdir"/final/"$f_name"nucalignCDStrans2.csv
  cat "$outdir"/final/*filteredfasta.csv | sed '1s Alignment,ID,error_type' >> "$outdir"/final/filteredfasta.csv
  virus=$(echo ""$f_name"")
  sixpack=$(grep -c "removed" "$wkd"/raw_seq/"$f_name"removedtoomanystops.csv)
  CDS=$(grep -c "removed" "$wkd"/raw_seq/"$f_name"removedCDS.csv)
  unknown=$(grep -c "removed" "$wkd"/raw_seq/"$f_name"removedmissingAA.csv)
  totalstart=$(grep -c ">" "$wkd"/input/"$f_name".fasta)
  totalused=$(wc -l "$file_in")
  percentseq=$(expr "$totalstart" / "$totalused" * "100")
  numofsites=$(wc -l "$file_out2")
  #variablesites=
  if [ ! -f "$outdir"/final/alignmentsummary.csv ];
  then
    echo "virus_name,sixpack_error,CDSnotfound,toomanyunknownAA,seqdownloaded,seqfiltered,seqused,percentused,totalsitesalignment,variable_sites" >> "$outdir"/final/alignmentsummary.csv
  fi
  echo ""$f_name","$sixpack","$CDS","$unknown","$totalstart","$totalfiltered","$totalused","$percentseq","$numofsites","$variablesites"" >> "$outdir"/final/alignmentsummary.csv

done
rm -rf "$outdir"/temp
