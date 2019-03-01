next_status="SETUP:";

#Setup commands.
yell() { printf "[%s\n" "${0##*/}] $*" | tr -s / >&2; }
die() { yell "$*"; exit 111; }
try() { (yell "Attempting: $*" && "$@" )  || die "Failed to complete: $*"; }
probe(){ command -v "$@" >/dev/null 2>&1 && yell "${next_status} Found $*" at \"$(which "$*")\" || die >&2 "${next_status} This script requires \""$*"\", but it's not installed."; }
realpath() { [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"; }
usage() { die "Usage: $0 [ -i INPUT (ab1, fasta, fastq) ]" 1>&2; }
exit_abnormal() { usage; exit 1; }

#Make sure all commands work
realpath "$0"
probe perl;
probe seqret;
probe blastn;
probe git;

body () {
  IFS= read -r header
  printf '%s\n' "$header"
  "$@"
}

OPTIONS=$(getopt -o :i:n:h -l input:,number:,help -- "$@")

usage() { die "Error: Usage: $0 -i [ab1/fastq/fasta] -n [max_target_seqs (>1)].";
	}

if [ $? -ne 0 ]; then
  	usage;
	exit 1
fi

eval set -- $OPTIONS

while true; do
  case "$1" in
    -i|--input) 	INPUT="$2";		shift ;;
    -n|--number)	max_target_seqs="$2"; 	shift ;;
    -h|--help)		usage ; 		shift ; exit 0 ;;
    --)       	        	 		shift ; break ;;
    *)        		usage ; 		exit 1 ;;
  esac
  shift
done

if [ ! -z "$@" ]; then 
	next_status="Warning :";

	yell "${next_status} \"$@\" is extraneous. Ignoring."
fi

if [ $INPUT != "ab1" ] && [ $INPUT != "fastq" ] && [ $INPUT != "fasta" ]; then

	next_status="Error:";
	die "${next_status} input flag [-i] must be ab1, fastq, or fasta.";

else

	yell "${next_status} input is set to ${INPUT}.";
fi

if [[ ! $max_target_seqs =~ ^[1-9][0-9]*$ ]]  ; then

	if [ -z $max_target_seqs ]; then 
		max_target_seqs="BLANK"; 
	fi
	
	next_status="Warning: "
	invalid_value=$max_target_seqs;
	max_target_seqs=10;

	yell "${next_status} Defaulting max_target_seqs to ${max_target_seqs}, was ${invalid_value}.";

else

	yell "${next_status} max_target_seqs is ${max_target_seqs}.";
	
fi

STRAIN_DEFINITIONS=/home/ryanward/Dropbox/Pietrasiak/JGI_strains.csv

trim=`ls -1 **/*trimmomatic*jar  2>/dev/null | wc -l | tr -d ' '`
if [ $trim = 0 ];
then
  yell "Trimmomatic not found...Getting it for you.";

  try wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip;
  try unzip Trimmomatic-0.38.zip;
  trimfile=$(realpath $(ls **/*trim*jar));
  yell "Trimmomatic is now at \"${trimfile}\"."
  try rm Trimmomatic*.zip;
else
  trimfile=$(realpath $(ls **/*trim*jar));
  yell "Found trimmomatic at \"${trimfile}\".";
fi
if [ $INPUT = "ab1" ] ; then
  step="EMBOSS:";
  count=`ls -1 *.ab1  2>/dev/null | wc -l | tr -d ' '`
  if [ $count = 0 ]

  then
    die "${step} No .ab1 file extensions in directory. Move them here to get started."
  else

    yell "${step} Found ${count} files with .ab1 extension."
    for x in $(ls *ab1);
    do
      next_status="Found:";
      outfile="$(basename "$x" .ab1).fq"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" = "0" ] ; then
          die "${step} ${x} file is zero length, incomplete or corrupt data.";
        else
          yell "${step} Converting ${x} to FASTQ format: ${outfile}" && seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} > ${outfile};
          next_status="Completed:";	
	fi
      fi
      outsize=$(cat $outfile | wc -l | tr -d ' ');
      if [ "$outsize" = "0" ] ; then
        yell "${step} ${next_status} ${outfile}: incomplete or corrupt... Re-trimming.";
        seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} > ${outfile};
      else
        yell "${step} ${next_status} ${outfile}, size: ${outsize} lines.";
	
      fi
    done;
  fi
fi

if [ $INPUT = "ab1" ] || [ $INPUT = "fastq" ] ; then
  step="TRIMMOMATIC:"
  count=`ls -1 *.fq 2>/dev/null | grep -v clean.fq  | wc -l | tr -d ' '`
  trimcount=0;
  if [ $count = 0 ]

  then
    die "${step} No .fq file extensions in directory."
  else

    yell "${step} Found ${count} files with .fq extension."
    yell "${step} Trimming ${count} files, ends with p-value of 0.05.";
    for x in $(ls *fq | grep -v clean.fq);

    do
      next_status="Found:";
      outfile="$(basename "$x" .fq).clean.fq"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" = "0" ] ; then
          die "${step} ${x}: corrupt input. Delete or move to continue.";
        else
      	  next_status="Completed:";
          java -jar ${trimfile} SE -phred64 ${x} ${outfile} LEADING:14 TRAILING:14 2>/dev/null;
        fi
      fi
      outsize=$(cat $outfile | wc -l | tr -d ' ');
      if [ "$outsize" = "0" ] ; then
        yell "${step} ${next_status} ${outfile}: incomplete or corrupt... Re-trimming.";
         java -jar ${trimfile} SE -phred64 ${x} ${outfile} LEADING:14 TRAILING:14 2>/dev/null;
        newoutsize=$(cat $outfile | wc -l | tr -d ' ');
      	next_status="Completed:";
        if [ "$newoutsize" = "0" ] ; then
          next_status="Warning:";  
          yell "${step} ${next_status} ${outfile} did not survive trimming."
        fi
      else
        yell "${step} ${next_status} ${outfile}, size: ${outsize} lines.";
          trim_count=$((trim_count+1));
      fi
    done;
  fi
  step="PASTE:"
  count=`ls -1 *.clean.fq  2>/dev/null | wc -l | tr -d ' '`
  if [ $count = 0 ]

  then
    die "${step} Found 0 trimmed reads. Aborting!"
  else

    yell "${step} Found ${trim_count} out of ${count}  that survived trimming, low quality files are zero-length."
    yell "${step} Converting from FASTQ to FASTA, all reads below 0.05 p-value."
    for x in $(ls *clean.fq);

    do
      outfile="$(basename "$x" .clean.fq).fa"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ $insize != "0" ];
      then
	      paste - - - - < ${x} | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${outfile};
      fi
    done;
  fi
fi
if [ $INPUT = "ab1" ] || [ $INPUT = "fastq" ]  || [ $INPUT = "fasta" ]; then
  step="BLASTN:";
  count=`ls -1 *.fa  2>/dev/null | wc -l | tr -d ' '`
  if [ $count = 0 ]

  then
    die "${step} No .fa file extensions in directory. Aborting!"
  else

    yell "${step} Found ${count} files with .fa extension."
    for x in $(ls *.fa);

    do
      outfile="$(basename "$x" .fa).blast.tsv"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" = "0" ] ; then
          yell "${step} ${x} Fasta size is zero.";
        else

          yell "${step} Blasting ${x}.";
          blastn -db nt -query $x -remote -max_target_seqs=${max_target_seqs} -out $outfile -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" || yell "${step} Timed out Blasting ${outfile}";
        fi
      else

        outsize=$(cat $outfile | wc -l | tr -d ' ');
        if [ "$outsize" = "0" ] || [ "$outsize" -lt $max_target_seqs ] ; then
          yell "${step} Found ${outfile}, size: ${outsize} lines. Attempting to fix."
          yell "${step} Blasting ${x}.";
          blastn -db nt -query $x -remote -max_target_seqs=${max_target_seqs} -out $outfile -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" || yell "${step} Timed out Blasting ${outfile}";
          newoutsize=$(cat $outfile | wc -l | tr -d ' ');
          if [ "$newoutsize" = "0" ] ; then
            yell "${step} ${outfile} was unable to complete."
          else

            yell "${step} Completed ${outfile}, size: ${newoutsize} lines.";
          fi
        else
          yell "${step} Found ${outfile}, size: ${outsize} lines.";
        fi
      fi
    done;
  fi
fi
  

  step="POST-PROCESS:"
  count=`ls -1 *.tsv 2>/dev/null | wc -l | tr -d ' '`
  fa_count=`ls -1 *.fa 2>/dev/null | wc -l | tr -d ' '`
  next_status="Found:"

  if [ $count != $fa_count ]

  then
    next_status="Error:";
    die "${step} ${next_status} Quantity mismatch of output (.tsv) and input (.fa) files. Please rerun this script."
  else
    for x in $(ls *tsv);
    do
      outfile="$(basename "$x" .tsv).fa"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" -lt "$max_target_seqs" ] ; then
          next_status="Error:";
	  die "${step} ${x} output (.tsv) is not ${max_target_seqs} in length, incomplete or corrupt data. Rerun this script.";
	fi
      fi
    done;
    
yell "${step} ${next_status} ${count} files with .tsv extension."
    echo "#q_sampleid,#q_primer,#q_seqid,#s_title,#s_acc,#s_seqid,#pident,#q_length,#s_length,#evalue,#bitscore,#q_genus,#q_species,#q_strain" > blast_out.csv;
    yell "${step} Concatenating results."
    try cat *tsv 2>/dev/null | sed 's/\,//g' | sed 's/\;//g' | sed 's/	/,/g' >> blast_out.csv;
    yell "${step} Blast results located at \"blast_out.csv\".";
  fi

###TODO###
#yell "${step} Extracting relevant information from results ..."
#try cat blast_out.csv 2>/dev/null | try body awk -vFS=, -vOFS=, '(NR!=1){match($1,/_([0-9]{1,2})[A-z]?_Pri/,sample);match ($1,/(Primer.*)/,primer); { print sample[1], primer[1],$0}}' > tmp.csv
#yell $"${step} Using ${STRAIN_DEFINITIONS} as source to extract query submission genus and species ..."
#try awk -vFS=, -vOFS=, '(NR==FNR){gen[$1]=$2;spec[$1]=$3;strain[$1]=$4;next;} (NR!=FNR) {if(FNR==1) print $0} {if(FNR!=1){print $0,gen[$1],spec[$1],strain[$1]} }' ${STRAIN_DEFINITIONS} tmp.csv > blast_out.csv 2>/dev/null && try rm tmp.csv 2>/dev/null
#paste -d, <( ${selefile} "col=q_sampleid,q_primer,q_genus,q_species,q_strain" blast_out.csv) <( ${selefile}  "col=s_" blast_out.csv) > tmp.csv
