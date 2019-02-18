#!/bin/bash


yell() { printf "[%s\n" "${0##*/}] $*" | tr -s / >&2; }
die() { yell "$*"; exit 111; }
try() { (yell "attempting: $*" && "$@" )  || die "failed to complete: $*"; }
probe(){ command -v "$@" >/dev/null 2>&1 && yell "Found $*" at \"$(which "$*")\" ...||  die >&2 "This script requires \""$*"\", but it's not installed. Aborting ..."; }

body () {
	IFS= read -r header
	printf '%s\n' "$header"
	"$@"
}


STRAIN_DEFINITIONS=/home/ryanward/Dropbox/Pietrasiak/JGI_strains.csv

probe seqret;

probe blastn;

probe git; 

sel=`ls -1 **/sel.awk  2>/dev/null | wc -l`

if [ $sel = 0 ];

then
	yell "Sel not found ... Getting it for you ...";

	try git clone https://github.com/ryandward/sel.git;

	selfile=$(realpath $(ls **/sel.awk));	

	yell "Sel is now at \"${selfile}\" ..."
else
	selfile=$(realpath $(ls **/sel.awk));	

	yell "Found sel at \"${selfile}\" ...";
fi;


trim=`ls -1 **/*trimmomatic*jar  2>/dev/null | wc -l`

if [ $trim = 0 ];

then
	yell "Trimmomatic not found ... Getting it for you ...";

	try wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip;

	try unzip Trimmomatic-0.38.zip;

	trimfile=$(realpath $(ls **/*trim*jar));	

	yell "Trimmomatic is now at \"${trimfile}\" ..."

	try rm Trimmomatic*.zip;
else
	trimfile=$(realpath $(ls **/*trim*jar));	

	yell "Found trimmomatic at \"${trimfile}\" ...";
fi;

step="EMBOSS:";

count=`ls -1 *.ab1  2>/dev/null | wc -l`

if [ $count = 0 ]

then 

	die "${step} No .ab1 file extensions in directory. Move them here to get started."
else 

	yell "${step} Found ${count} files with .ab1 extension ..."

	for x in $(ls *ab1); 

	do 
		outfile=$x.fq;

		insize=$(stat -c%s $x);

		if [ ! -f $outfile ] ; then

			if [ "$insize" = "0" ] ; then

				die "${step} ${x} file is zero length, incomplete or corrupt data.";

			else 		
				yell "${step} Converting ${x} to FASTQ format: ${outfile}" && try seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} > ${x}.fq;
			fi


		fi

		outsize=$(stat -c%s $outfile);

		if [ "$outsize" = "0" ] ; then

			yell "${step} ${outfile} is zero length, incomplete or corrupt conversion. Re-attempting ...";

			try seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} > ${x}.fq;

		else
			yell "${step} Found ${outfile} size: ${outsize} ...";
		fi;

	done;
fi;


step="TRIMMOMATIC:"

count=`ls -1 *.ab1.fq  2>/dev/null | wc -l`

if [ $count = 0 ]

then 

	die "${step} No .ab1.fq file extensions in directory."
else 

	yell "${step} Found ${count} files with .ab1.fq extension ..."

	yell "${step} Trimming ${count} files, ends with p-value of 0.05 ...";

	for x in $(ls *ab1.fq); 

	do 
		outfile="$(basename "$x" .ab1.fq).clean.fq"

		insize=$(stat -c%s $x);

		if [ ! -f $outfile ] ; then

			if [ "$insize" = "0" ] ; then

				die "${step} ${x} file is zero length, incomplete or corrupt data.";

			else 		
				try java -jar ${trimfile} SE -phred64 ${x} ${outfile} LEADING:14 TRAILING:14 2>/dev/null;

			fi


		fi

		outsize=$(stat -c%s $outfile);

		if [ "$outsize" = "0" ] ; then

			yell "${step} ${outfile} is zero length, incomplete, corrupt conversion, or very low quality data ... Re-attempting trimmomatic ...";

			try java -jar ${trimfile} SE -phred64 ${x} ${outfile} LEADING:14 TRAILING:14 2>/dev/null;


			newoutsize=$(stat -c%s $outfile);
			if [ "$newoutsize" = "0" ] ; then

				yell "${step} ${outfile} did not survive trimming ..."


				count=$((count-1));

			fi
		else
			yell "${step} Found ${outfile} size: ${outsize} ...";
		fi;

	done;
fi;


pretrimcount=$count;


step="PASTE:"

count=`ls -1 *.clean.fq  2>/dev/null | wc -l`

if [ $count = 0 ]

then 

	die "${step} Found 0 trimmed reads out of ${pretrimcount} ... Aborting!"
else 

	yell "${step} Found ${pretrimcount} out of ${count}  that survived trimming, low quality files are zero-length ..."


	yell "${step} Converting from FASTQ to FASTA, all reads above 0.05 p-value."

	for x in $(ls *clean.fq); 


	do 
		outfile="$(basename "$x" .fq).fa"
		insize=$(stat -c%s $x);

		if [ $insize != "0" ];

		then

			paste - - - - < ${x} | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${outfile}; 

		fi
	done;


fi

step="BLASTN:";

count=`ls -1 *.clean.fa  2>/dev/null | wc -l`

if [ $count = 0 ]

then 

	die "${step} No .clean.fa file extensions in directory. Something went wrong ... Aborting!"
else 

	yell "${step} Found ${count} files with .clean.fa extension ..."

	for x in $(ls *clean.fa);


	do

		outfile="$(basename "$x" .clean.fa).blast.tsv"

		insize=$(stat -c%s $x);

		if [ ! -f $outfile ] ; then

			if [ "$insize" = "0" ] ; then


				yell "${step} ${x} FASTA size zero --> QC fail.";

			else 		

				yell "${step} Blasting ${x} ...";

				timeout --foreground 5m blastn -db nt -query $x -remote -max_target_seqs=10 -out $outfile -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" || yell "${step} Timed out Blasting ${outfile}";

			#	yell "${step} Exit code: $?"; 
			
			fi


		else	

			outsize=$(stat -c%s $outfile);

			if [ "$outsize" = "0" ] ; then
				yell "${step} Found ${outfile} size: ${outsize}. Attempting to fix ... "
		
				yell "${step} Blasting ${x} ...";

				timeout --foreground 5m blastn -db nt -query $x -remote -max_target_seqs=10 -out $outfile -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" || yell "${step} Timed out Blasting ${outfile}";

			#	yell "${step} Exit code: $?"; 
			
				newoutsize=$(stat -c%s $outfile);
			
				if [ "$newoutsize" = "0" ] ; then

				
					yell "${step} ${outfile} was unable to complete ... "
				else
				
					yell "${step} Completed ${outfile} size: ${newoutsize} ...";

				fi;
			else
				yell "${step} Found ${outfile} size: ${outsize} ...";
			fi;
		fi;

	done;
fi;


step="POST-PROCESS:"

echo "#q_sampleid,#q_primer,#q_seqid,#s_title,#s_acc,#s_seqid,#pident,#q_length,#s_length,#evalue,#bitscore,#q_genus,#q_species,#q_strain" > blast_out.csv;

yell "${step} Concatenating results ..."
try cat *tsv 2>/dev/null | sed 's/\,//g' | sed 's/\;//g' | sed 's/\t/,/g' >> blast_out.csv;

yell "${step} Extracting relevant information from results ..."
try cat blast_out.csv 2>/dev/null | try body awk -vFS=, -vOFS=, '(NR!=1){match($1,/_([0-9]{1,2})[A-z]?_Pri/,sample);match ($1,/(Primer.*)/,primer); { print sample[1], primer[1],$0}}' > tmp.csv 


yell $"${step} Using ${STRAIN_DEFINITIONS} as source to extract query submission genus and species ..."
try awk -vFS=, -vOFS=, '(NR==FNR){gen[$1]=$2;spec[$1]=$3;strain[$1]=$4;next;} (NR!=FNR) {if(FNR==1) print $0} {if(FNR!=1){print $0,gen[$1],spec[$1],strain[$1]} }' ${STRAIN_DEFINITIONS} tmp.csv > blast_out.csv 2>/dev/null && try rm tmp.csv 2>/dev/null

paste -d, <( ${selfile} "col=q_sampleid,q_primer,q_genus,q_species,q_strain" blast_out.csv) <( ${selfile}  "col=s_" blast_out.csv) > tmp.csv 

