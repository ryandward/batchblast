step="Setup:";
#Setup commands.
yell() { printf ["%s\n" "${0##*/}] ${step} ${status} $*" | tr -s / >&2; }
die() { yell "$*"; exit 111; }
try() { (yell "Trying: $*" && "$@" )  || die "Failed: $*"; }
probe(){ command -v "$@" >/dev/null 2>&1 && yell "Found: $*" || die >&2 "This script requires \""$*"\", but it's not installed."; }
realpath() { [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"; }
usage() { die "Usage: $0 [ -i INPUT (ab1, fasta, fastq) ]" 1>&2; }
exit_abnormal() { usage; exit 1; }

#Make sure all commands work
status="Probe:";
realpath "$0"
probe perl;
probe seqret;
probe blastn;
probe git;
probe timeout;

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
	status="Warning :";
	yell "\"$@\" is extraneous. Ignoring."
fi

if [ -z $INPUT ]; then
	status="Error :";
	die "input flag [-i] cannot be blank, and must be ab1, fastq, or fasta.";
fi

if [ $INPUT != "ab1" ] && [ $INPUT != "fastq" ] && [ $INPUT != "fasta" ]; then
	status="Error:";
	die "input flag [-i] must be ab1, fastq, or fasta.";

else

	status="Variable:";
	yell "input set to ${INPUT}.";
fi

if [[ $max_target_seqs -gt 100 ]]; then
	status="Warning:";
	yell $max_target_seqs;
	max_target_seqs=100
	yell "maximum returned results must be between 1 and ${max_target_seqs}. Setting to ${max_target_seqs}.";
fi

if [[ ! $max_target_seqs =~ ^[1-9][0-9]*$ ]]  ; then

	if [ -z $max_target_seqs ]; then 
		max_target_seqs="BLANK"; 
	fi
	
	status="Warning:"
	invalid_value=$max_target_seqs;
	max_target_seqs=10;

	yell "defaulting max_target_seqs to ${max_target_seqs}, was ${invalid_value}.";

else

	yell "max_target_seqs is ${max_target_seqs}.";
	
fi

STRAIN_DEFINITIONS=/home/ryanward/Dropbox/Pietrasiak/JGI_strains.csv

trim=`ls -1 **/*trimmomatic*jar  2>/dev/null | wc -l | tr -d ' '`
if [ $trim = 0 ];
then
  status="Warning:";
  yell "Trimmomatic not found... Installing.";
  try wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip;
  try unzip Trimmomatic-0.38.zip;
  trimfile=$(realpath $(ls **/*trim*jar));
  Status="Success:";
  yell "Trimmomatic is now at \"${trimfile}\"."
  try rm Trimmomatic*.zip;
else
  status="Found";
  trimfile=$(realpath $(ls **/*trim*jar));
  yell "Trimmomatic at \"${trimfile}\".";
fi
if [ $INPUT = "ab1" ] ; then
  step="EMBOSS:";
  count=`ls -1 *.ab1  2>/dev/null | wc -l | tr -d ' '`
  if [ $count = 0 ]

  then
    status="Error:";
    die "No .ab1 file extensions in directory. Move them here to get started."
  else
    status="Found:":
    yell "${count} files with .ab1 extension."
    for x in $(ls *ab1);
    do
      status="Found:";
      outfile="$(basename "$x" .ab1).fq"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" = "0" ] ; then
          status="Error:";
          die "${x} file is zero length, incomplete or corrupt data.";
        else
          status="Converting:";
          yell "${x} to FastQ format: ${outfile}" && seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} > ${outfile};
          status="Completed:";	
	fi
      fi
      outsize=$(cat $outfile | wc -l | tr -d ' ');
      if [ "$outsize" = "0" ] ; then
        yell "${outfile}: incomplete or corrupt... Re-trimming.";
        seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} > ${outfile};
      else
        yell "${outfile}, size: ${outsize} lines.";
	
      fi
    done;
  fi
fi

if [ $INPUT = "ab1" ] || [ $INPUT = "fastq" ] ; then
  step="Trimmomatic:"
  status="Begin:"
  count=`ls -1 *.fq 2>/dev/null | grep -v trim.fq  | wc -l | tr -d ' '`
  trimcount=0;
  if [ $count = 0 ]; then
    status="Error:";
    die "No .fq file extensions in directory."
  else
    status="Found:";
    yell "${count} files with .fq extension."
    status="Trimming:";
    yell "${count} files, ends with p-value of 0.05.";
    for x in $(ls *fq | grep -v trim.fq);

    do
      status="Found:";
      outfile="$(basename "$x" .fq).trim.fq"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" = "0" ] ; then
          die "${x}: corrupt input. Delete or move to continue.";
        else
      	  status="Completed:";
          java -jar ${trimfile} SE -phred64 ${x} ${outfile} LEADING:14 TRAILING:14 2>/dev/null;
        fi
      fi
      outsize=$(cat $outfile | wc -l | tr -d ' ');
      if [ "$outsize" = "0" ] ; then
        yell "${outfile}: incomplete or corrupt... Re-trimming.";
         java -jar ${trimfile} SE -phred64 ${x} ${outfile} LEADING:14 TRAILING:14 2>/dev/null;
        newoutsize=$(cat $outfile | wc -l | tr -d ' ');
      	status="Completed:";
        if [ "$newoutsize" = "0" ] ; then
          status="Warning:";  
          yell "${outfile} did not survive trimming."
        fi
      else
        yell "${outfile}, size: ${outsize} lines.";
          trim_count=$((trim_count+1));
      fi
    done;
  fi
  step="Data Mining:"
  count=`ls -1 *.trim.fq  2>/dev/null | wc -l | tr -d ' '`
  if [ $count = 0 ]

  then
    status="Error:";
    die "0 trimmed reads. Aborting!"
  else
    status="Found:";
    yell "${trim_count} out of ${count} that survived trimming."
    status="Converting:";
    yell "FastQ to Fasta, p-value < 0.05.";
    for x in $(ls *trim.fq);

    do
      outfile="$(basename "$x" .trim.fq).fa"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ $insize != "0" ];
      then
        paste - - - - < ${x} | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${outfile};
      fi
    done;
  fi
fi
if [ $INPUT = "ab1" ] || [ $INPUT = "fastq" ]  || [ $INPUT = "fasta" ]; then
  step="BlastN:";
  count=`ls -1 *.fa  2>/dev/null | wc -l | tr -d ' '`
  if [ $count = 0 ]

  then
    status="Error:";
    die "No .fa file extensions in directory. Aborting!"
  else
    status="Found:";
    yell "${count} files with .fa extension."
    for x in $(ls *.fa);

    do
      outfile="$(basename "$x" .fa).blast.tsv"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" = "0" ] ; then
          yell "${x} Fasta size is zero.";
        else
          status="Blasting:";
          yell "${x}.";
          timeout --foreground 2m blastn -db nt -query $x -remote -max_target_seqs=${max_target_seqs} -out $outfile -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" || (status="Warning:" && yell "Operation timed out.";)
        fi
      else

        outsize=$(cat $outfile | wc -l | tr -d ' ');
        if [ "$outsize" = "0" ] || [ "$outsize" -lt $max_target_seqs ] ; then
          status="Found:";
          yell "${outfile} size: ${outsize} lines. Attempting to fix."
          
          status="Blasting:";
          yell "${x}.";
          timeout --foreground 2m blastn -db nt -query $x -remote -max_target_seqs=${max_target_seqs} -out $outfile -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" || (status="Warning:" && yell "Operation timed out.";)
          newoutsize=$(cat $outfile | wc -l | tr -d ' ');
          if [ "$newoutsize" = "0" ] ; then
            yell "${outfile} was unable to complete."
          else

            yell "Completed ${outfile}, size: ${newoutsize} lines.";
          fi
        else
          yell "Found ${outfile}, size: ${outsize} lines.";
        fi
      fi
    done;
  fi
fi
  

  step="POST-PROCESS:"
  count=`ls -1 *.tsv 2>/dev/null | wc -l | tr -d ' '`
  fa_count=`ls -1 *.fa 2>/dev/null | wc -l | tr -d ' '`
  status="Found:"

  if [ $count != $fa_count ]

  then
    status="Error:";
    die "Quantity mismatch of output (.tsv) and input (.fa) files. Please rerun this script."
  else
    for x in $(ls *tsv);
    do
      outfile="$(basename "$x" .tsv).fa"
      insize=$(cat $x | wc -l | tr -d ' ');
      if [ ! -f $outfile ] ; then
        if [ "$insize" -lt "$max_target_seqs" ] ; then
          status="Error:";
	  die "${x} output (.tsv) is not ${max_target_seqs} in length, incomplete or corrupt data. Rerun this script.";
	fi
      fi
    done;
    
yell "${count} files with .tsv extension."
    echo "#q_sampleid,#q_primer,#q_seqid,#s_title,#s_acc,#s_seqid,#pident,#q_length,#s_length,#evalue,#bitscore,#q_genus,#q_species,#q_strain" > blast_out.csv;
    yell "Concatenating results."
    try cat *tsv 2>/dev/null | sed 's/\,//g' | sed 's/\;//g' | sed 's/	/,/g' >> blast_out.csv;
    yell "Blast results located at \"blast_out.csv\".";
  fi

###TODO###
#yell "Extracting relevant information from results ..."
#try cat blast_out.csv 2>/dev/null | try body awk -vFS=, -vOFS=, '(NR!=1){match($1,/_([0-9]{1,2})[A-z]?_Pri/,sample);match ($1,/(Primer.*)/,primer); { print sample[1], primer[1],$0}}' > tmp.csv
#yell $"Using ${STRAIN_DEFINITIONS} as source to extract query submission genus and species ..."
#try awk -vFS=, -vOFS=, '(NR==FNR){gen[$1]=$2;spec[$1]=$3;strain[$1]=$4;next;} (NR!=FNR) {if(FNR==1) print $0} {if(FNR!=1){print $0,gen[$1],spec[$1],strain[$1]} }' ${STRAIN_DEFINITIONS} tmp.csv > blast_out.csv 2>/dev/null && try rm tmp.csv 2>/dev/null
#paste -d, <( ${selefile} "col=q_sampleid,q_primer,q_genus,q_species,q_strain" blast_out.csv) <( ${selefile}  "col=s_" blast_out.csv) > tmp.csv
