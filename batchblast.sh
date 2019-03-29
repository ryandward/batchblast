WORKDIR="/home/ryanward/Dropbox/Hutten/Blast_Fun"
STEP="Setup:"
#Setup commands.
yell() { printf ["%s\n" "${0##*/}] ${STEP} ${STATUS} $*" | tr -s / >&2; }
die() {
  yell "$*"
  exit 111
}
try() { (yell "Trying: $*" && "$@") || die "Failed: $*"; }
probe() { command -v "$@" >/dev/null 2>&1 && yell "Found: $*" || die >&2 "This script requires \""$*"\", but it's not installed."; }
realpath() { [[ $1 == /* ]] && echo "$1" || echo "$PWD/${1#./}"; }
usage() { die "Usage: $0 [ -i INPUT (ab1, fasta, fastq) ]" 1>&2; }
exit_abnormal() {
  usage
  exit 1
}

#Make sure all commands work
STATUS="Probe:"
realpath "$0"
probe perl
probe seqret
probe blastn
probe git
probe timeout

body() {
  IFS= read -r header
  printf '%s\n' "$header"
  "$@"
}

OPTIONS=$(getopt -o :i:n:h -l input:,number:,help -- "$@")

usage() {
  die "Error: Usage: $0 -i [ab1/fastq/fasta] -n [MAX_TARGET_SEQS (>1)]."
}

if [ $? -ne 0 ]; then
  usage
  exit 1
fi

eval set -- $OPTIONS

while true; do
  case "$1" in
  -i | --input)
    INPUT="$2"
    shift
    ;;
  -n | --number)
    MAX_TARGET_SEQS="$2"
    shift
    ;;
  -h | --help)
    usage
    shift
    exit 0
    ;;
  --)
    shift
    break
    ;;
  *)
    usage
    exit 1
    ;;
  esac
  shift
done

if [ ! -z "$@" ]; then
  STATUS="Warning :"
  yell "\"$@\" is extraneous. Ignoring."
fi

if [ -z $INPUT ]; then
  STATUS="Error :"
  die "input flag [-i] cannot be blank, and must be ab1, fastq, or fasta."
fi

if [ $INPUT != "ab1" ] && [ $INPUT != "fastq" ] && [ $INPUT != "fasta" ]; then
  STATUS="Error:"
  die "input flag [-i] must be ab1, fastq, or fasta."

else

  STATUS="Variable:"
  yell "input set to ${INPUT}."
fi

if [[ $MAX_TARGET_SEQS -gt 100 ]]; then
  STATUS="Warning:"
  yell $MAX_TARGET_SEQS
  MAX_TARGET_SEQS=100
  yell "maximum returned results must be between 1 and ${MAX_TARGET_SEQS}. Setting to ${MAX_TARGET_SEQS}."
fi

if [[ ! $MAX_TARGET_SEQS =~ ^[1-9][0-9]*$ ]]; then

  if [ -z $MAX_TARGET_SEQS ]; then
    MAX_TARGET_SEQS="BLANK"

  fi

  STATUS="Warning:"
  invalid_value=$MAX_TARGET_SEQS
  MAX_TARGET_SEQS=10

  yell "defaulting MAX_TARGET_SEQS to ${MAX_TARGET_SEQS}, was ${invalid_value}."

else
  yell "MAX_TARGET_SEQS is ${MAX_TARGET_SEQS}."

fi

STRAIN_DEFINITIONS=/home/ryanward/Dropbox/Pietrasiak/JGI_strains.csv

trim=$(ls -1 **/*trimmomatic*jar 2>/dev/null | wc -l | tr -d ' ')

if [ $trim = 0 ]; then
  STATUS="Warning:"
  yell "Trimmomatic not found... Installing."
  try wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
  try unzip Trimmomatic-0.38.zip
  TRIMFILE=$(realpath $(ls **/*trim*jar))
  STATUS="Success:"
  yell "Trimmomatic is now at \"${TRIMFILE}\"."
  try rm Trimmomatic*.zip

else
  STATUS="Found"
  TRIMFILE=$(realpath $(ls **/*trim*jar))
  yell "Trimmomatic at \"${TRIMFILE}\"."

fi

cd "$WORKDIR";

for x in $(ls); do dos2unix $x; done;

if [ $INPUT = "ab1" ]; then
  STEP="EMBOSS:"
  PROGRESS=0;
  COUNT=$(ls -1 *.ab1 2>/dev/null | wc -l | tr -d ' ')

  if [ $COUNT = 0 ]; then
    STATUS="Error:"
    die "No .ab1 file extensions in directory. Move them here to get started."

  else
    STATUS="Found:":
    yell "${COUNT} files with .ab1 extension."

    for x in $(ls *ab1); do
      ((PROGRESS++))
      STATUS="($PROGRESS/$COUNT) Found:"
      OUTFILE="$(basename "$x" .ab1).fq"
      INSIZE=$(cat $x | wc -l | tr -d ' ')

      if [ ! -f $OUTFILE ]; then

        if [ "$INSIZE" = "0" ]; then
          STATUS="($PROGRESS/$COUNT) Warning:"
          die "${x} file is zero length, incomplete or corrupt data."

        else
          STATUS="($PROGRESS/$COUNT) Converting:"
          yell "${x} to FastQ format: ${OUTFILE}" && seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} >${OUTFILE}
          STATUS="($PROGRESS/$COUNT) Completed:"

        fi

      fi

      OUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')

      if [ "$OUTSIZE" = "0" ]; then
        yell "${OUTFILE}: incomplete or corrupt... Re-trimming."
        seqret -sformat abi -osformat fastq -auto -stdout -sequence ${x} >${OUTFILE}

      else
        yell "${OUTFILE}, size: ${OUTSIZE} lines."

      fi

    done

  fi

fi

if [ $INPUT = "ab1" ] || [ $INPUT = "fastq" ]; then
  STEP="Trimmomatic:"
  PROGRESS=0;
  STATUS="Begin:"
  COUNT=$(ls -1 *.fq 2>/dev/null | grep -v trim.fq | wc -l | tr -d ' ')
  trimCOUNT=0

  if [ $COUNT = 0 ]; then
    STATUS="Error:"
    die "No .fq file extensions in directory."

  else
    STATUS="Found:"
    yell "${COUNT} files with .fq extension."
    STATUS="Initialize:"
    yell "${COUNT} files, ends with p-value of 0.05."

    for x in $(ls *fq | grep -v trim.fq); do
      ((PROGRESS++))
      STATUS="($PROGRESS/$COUNT) Found:"
      OUTFILE="$(basename "$x" .fq).trim.fq"
      INSIZE=$(cat $x | wc -l | tr -d ' ')

      if [ ! -f $OUTFILE ]; then

        if [ "$INSIZE" = "0" ]; then
          die "${x}: corrupt input. Delete or move to continue."

        else
          STATUS="($PROGRESS/$COUNT) Completed:"
          java -jar ${TRIMFILE} SE -phred64 ${x} ${OUTFILE} LEADING:14 TRAILING:14 2>/dev/null

        fi

      fi

      OUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')

      if [ "$OUTSIZE" = "0" ]; then
        yell "${OUTFILE}: incomplete or corrupt... Re-trimming."
        java -jar ${TRIMFILE} SE -phred64 ${x} ${OUTFILE} LEADING:14 TRAILING:14 2>/dev/null
        newOUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')
        STATUS="($PROGRESS/$COUNT) Completed:"

        if [ "$newOUTSIZE" = "0" ]; then
          STATUS="($PROGRESS/$COUNT) Warning:"
          yell "${OUTFILE} did not survive trimming."

        fi

      else
        yell "${OUTFILE}, size: ${OUTSIZE} lines."
        trim_COUNT=$((trim_COUNT + 1))

      fi

    done

  fi

  STEP="Formatting:"
  PROGRESS=0;
  COUNT=$(ls -1 *.trim.fq 2>/dev/null | wc -l | tr -d ' ')

  if [ $COUNT = 0 ]; then
    STATUS="Error:"
    die "0 trimmed reads. Aborting!"

  else
    STATUS="Found:"

    yell "${trim_COUNT} out of ${COUNT} that survived trimming."
    COUNT=$trim_COUNT;

    yell "FastQ to FASTA for p-value < 0.05."

    for x in $(ls *trim.fq); do
      ((PROGRESS++))
      STATUS="($PROGRESS/$COUNT) Converting:"

      OUTFILE="$(basename "$x" .trim.fq).fa"
      INSIZE=$(cat $x | wc -l | tr -d ' ')
      yell "${x} --> ${OUTFILE}"

      if [ $INSIZE != "0" ]; then
        paste - - - - <${x} | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" >${OUTFILE}

      fi

    done

  fi

fi

if [ $INPUT = "ab1" ] || [ $INPUT = "fastq" ] || [ $INPUT = "fasta" ]; then

  STEP="BlastN:"
  PROGRESS=0;

  COUNT=$(ls -1 *.fa 2>/dev/null | wc -l | tr -d ' ')

  if [ $COUNT = 0 ]; then
    STATUS="Error:"
    die "No .fa file extensions in directory. Aborting!"

  else
    STATUS="Found:"
    yell "${COUNT} files with .fa extension."

    for x in $(ls *.fa); do
      ((PROGRESS++))
      OUTFILE="$(basename "$x" .fa).blast.tsv"

      INSIZE=$(cat $x | wc -l | tr -d ' ')

      if [ ! -f $OUTFILE ]; then

        if [ "$INSIZE" = "0" ]; then
          yell "${x} Fasta size is zero."

        else
          STATUS="($PROGRESS/$COUNT) Blasting:"
          yell "${x}."
          timeout --foreground 2m blastn -db nt -query $x -remote -max_target_seqs=${MAX_TARGET_SEQS} -out $OUTFILE -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" 2>&1 |

            while read line; do
              STATUS="NCBI returned:"
              yell $line

            done

            STATUS="($PROGRESS/$COUNT) Blasting:"

        fi

      else

        OUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')

        if [ "$OUTSIZE" = "0" ] || [ "$OUTSIZE" -lt $MAX_TARGET_SEQS ]; then
          STATUS="($PROGRESS/$COUNT) Found:"
          yell "${OUTFILE} size: ${OUTSIZE} lines. Attempting to fix."
          STATUS="($PROGRESS/$COUNT) Blasting:"
          yell "${x}."
          timeout --foreground 2m blastn -db nt -query $x -remote -max_target_seqs=${MAX_TARGET_SEQS} -out $OUTFILE -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" 2>&1 |

            while read line; do
              STATUS="($PROGRESS/$COUNT) NCBI returned:"
              yell $line

            done

          STATUS="($PROGRESS/$COUNT) Blasting:"

          newOUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')

          if [ "$newOUTSIZE" = "0" ]; then
            STATUS="($PROGRESS/$COUNT) Warning:"
            yell "Unable to initialize ${OUTFILE}, probable NCBI network throttle."
            STATUS="($PROGRESS/$COUNT) Blasting:"

          else

            yell "Completed ${OUTFILE}, size: ${newOUTSIZE} lines."

          fi

        else
          STATUS="($PROGRESS/$COUNT) Found:"
          yell "${OUTFILE}, size: ${OUTSIZE} lines."
          STATUS="($PROGRESS/$COUNT) Blasting:"

        fi

      fi

    done

  fi

fi

STEP="POST-PROCESS:"
COUNT=$(ls -1 *.tsv 2>/dev/null | wc -l | tr -d ' ')
fa_COUNT=$(ls -1 *.fa 2>/dev/null | wc -l | tr -d ' ')
STATUS="Found:"

if [ $COUNT != $fa_COUNT ]; then
  STATUS="Error:"
  die "Quantity mismatch of output (.tsv) and input (.fa) files. Please rerun this script."

else

  for x in $(ls *tsv); do
    OUTFILE="$(basename "$x" .tsv).fa"
    INSIZE=$(cat $x | wc -l | tr -d ' ')

    if [ ! -f $OUTFILE ]; then

      if [ "$INSIZE" -lt "$MAX_TARGET_SEQS" ]; then

        STATUS="Error:"
        die "${x} output .tsv file is not ${MAX_TARGET_SEQS} in length, incomplete or corrupt data. Rerun this script."

      fi

    fi

  done

  for x in $(ls); do dos2unix $x; done;

  yell "${COUNT} files with .tsv extension."
  yell "Concatenating results."
  echo "Query,Subject Title,Subject Accession,Accession URL,FASTA URL,Percent Identical,Query Length,Subject Length,E Value,Bitscore" >blast_out.csv
  try cat *tsv 2>/dev/null | awk 'BEGIN { FS="\t"; OFS="," } {rebuilt=0; for(i=1; i<=NF; ++i) {if ($i ~ /,/ && $i !~ /^".*"$/) { $i = "\"" $i "\"";rebuilt=1 }}  if (!rebuilt) { $1=$1 }print}' | awk 'BEGIN{OFS=",";FPAT = "([^,]+)|(\"[^\"]+\")"}NR!=1{split ($4,pieces,"|"); $4="=HYPERLINK(\"https://www.ncbi.nlm.nih.gov/nuccore/"pieces[2]"\""",=HYPERLINK(\"https://www.ncbi.nlm.nih.gov/nuccore/"pieces[4]"?report=fasta\")"; $1=$1; print}' >> blast_out.csv
  yell "Blast results located at \"blast_out.csv\"."

fi

###TODO###
#yell "Extracting relevant information from results ..."
#try cat blast_out.csv 2>/dev/null | try body awk -vFS=, -vOFS=, '(NR!=1){match($1,/_([0-9]{1,2})[A-z]?_Pri/,sample);match ($1,/(Primer.*)/,primer); { print sample[1], primer[1],$0}}' > tmp.csv
#yell $"Using ${STRAIN_DEFINITIONS} as source to extract query submission genus and species ..."
#try awk -vFS=, -vOFS=, '(NR==FNR){gen[$1]=$2;spec[$1]=$3;strain[$1]=$4;next;} (NR!=FNR) {if(FNR==1) print $0} {if(FNR!=1){print $0,gen[$1],spec[$1],strain[$1]} }' ${STRAIN_DEFINITIONS} tmp.csv > blast_out.csv 2>/dev/null && try rm tmp.csv 2>/dev/null
#paste -d, <( ${selefile} "col=q_sampleid,q_primer,q_genus,q_species,q_strain" blast_out.csv) <( ${selefile}  "col=s_" blast_out.csv) > tmp.csv
