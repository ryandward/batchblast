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
BLASTPATH="$(realpath "$0")"
BLASTDIR="$(dirname "$BLASTPATH")"
ARGS="$@"

probe perl
probe seqret
probe blastn
probe git
probe timeout
probe dos2unix

body() {
  IFS= read -r header
  printf '%s\n' "$header"
  "$@"
}

OPTIONS=$(getopt -o :i:n:w:h -l input:,number:,workdir:,help -- "$@")

usage() {
  die "Error: Usage: $0 -i [ab1/fastq/fasta] -n [target sequences] -w [work directory]."
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
    -w | --workdir)
    WORKDIR="$2"
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

if [ -z $WORKDIR ]; then

  STATUS="Error:"
  die "Work directory flag [-w] cannot be blank."

else

  STATUS="Variable:"
  yell "Workdir set to ${WORKDIR}."
  ls "$WORKDIR" >/dev/null || die "Could not find $WORKDIR"
  for x in $(ls $WORKDIR); do dos2unix &>/dev/null $x; done;

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

trim=$(ls -1 $BLASTDIR/**/*trimmomatic*jar 2>/dev/null | wc -l | tr -d ' ')

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


if [ $INPUT = "ab1" ]; then
  STEP="EMBOSS:"
  PROGRESS=0;
  COUNT=$(ls -1 *.ab1 2>/dev/null | wc -l | tr -d ' ')

  if [ $COUNT = 0 ]; then
    STATUS="Error:"
    die "No .ab1 file extensions in directory. Move them here to get started."

  else
    STATUS="Found":
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
          timeout --foreground 3m blastn -db nt -query $x -remote -max_target_seqs=${MAX_TARGET_SEQS} -out $OUTFILE -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" 2>&1 |

          while read line; do
            STATUS="NCBI returned:"
            yell $line

            if [[ $line == *"failed"* ]] || [[ $line == *"Failed"* ]]; then
              ps -ef | grep $x | grep $USER | grep timeout | grep -v grep | awk '{print $2}' | xargs kill && break
            fi


          done

          yell "${OUTFILE}, size: ${newOUTSIZE} lines."
          yell "Pausing..."
          sleep 60

          STATUS="($PROGRESS/$COUNT) Blasting:"

        fi

      else

        OUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')

        if [ "$OUTSIZE" = "0" ] || [ "$OUTSIZE" -lt $MAX_TARGET_SEQS ]; then
          STATUS="($PROGRESS/$COUNT) Found:"
          yell "${OUTFILE} size: ${OUTSIZE} lines. Attempting to fix."
          STATUS="($PROGRESS/$COUNT) Blasting:"
          yell "${x}."
          timeout --foreground 3m blastn -db nt -query $x -remote -max_target_seqs=${MAX_TARGET_SEQS} -out $OUTFILE -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" 2>&1 |

          while read line; do
            STATUS="($PROGRESS/$COUNT) NCBI returned:"
            yell $line

            if [[ $line == *"failed"* ]] || [[ $line == *"Failed"* ]]; then
              ps -ef | grep $x | grep $USER | grep timeout | grep -v grep | awk '{print $2}' | xargs kill && break
            fi

          done
          yell "${OUTFILE}, size: ${newOUTSIZE} lines."
          yell "Pausing..."
          sleep 60

          STATUS="($PROGRESS/$COUNT) Blasting:"

          newOUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')

          if [ "$newOUTSIZE" = "0" ]; then
            STATUS="($PROGRESS/$COUNT) Warning:"
            yell "Unable to initialize ${OUTFILE}, probable NCBI network throttle."
            STATUS="($PROGRESS/$COUNT) Blasting:"

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

STEP="Data Query:"
COUNT=$(ls -1 *.tsv 2>/dev/null | wc -l | tr -d ' ')
fa_COUNT=$(ls -1 *.fa 2>/dev/null | wc -l | tr -d ' ')
STATUS="Found:"

if [ $COUNT != $fa_COUNT ]; then
  STATUS="Error:"
  exec $BLASTPATH $ARGS && die "Quantity mismatch of output (.tsv) and input (.fa) files. This script will rerun."

else

  for x in $(ls *tsv); do
    OUTFILE="$(basename "$x" .tsv).fa"
    INSIZE=$(cat $x | wc -l | tr -d ' ')

    if [ ! -f $OUTFILE ]; then

      if [ "$INSIZE" -lt "$MAX_TARGET_SEQS" ]; then

        STATUS="Error:"
        exec $BLASTPATH $ARGS && die "${x} output .tsv is smaller than expected. This script will rerun."

      fi

    fi

  done

  for x in $(ls); do dos2unix &>/dev/null $x; done;

  yell "${COUNT} files with .tsv extension."
  yell "Gathering data from taxonomy.jgi-psf.org, this may take a while."
  echo "Query Name","Target Accession","Target Genus","Target Species","Target Strain","Identical","Query Length","Target Length","E-Value","Bitscore","Fasta","Target Description","Reported Taxonomic Name" > blast_out.csv
  try cat *tsv 2>/dev/null |
  exec "$BLASTDIR/rebuild.awk" |
  exec "$BLASTDIR/extract.awk" |
  tee -a blast_out.csv
  yell "Blast results located at \"blast_out.csv\"."

fi
