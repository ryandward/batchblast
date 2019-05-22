#Setup commands.

STEP="Setup:"

yell() { printf ["%s\n" "${0##*/}] ${STEP} ${STATUS} $*" | tr -s / >&2; }
die() {
  yell "$*"
  exit 111
}
try() { (yell "Trying: $*" && "$@") || die "Failed: $*"; }
probe() { command -v "$@" >/dev/null 2>&1 && yell "Found: $*" || die >&2 "This script requires \""$*"\", but it's not installed."; }
realpath() { [[ $1 == /* ]] && echo "$1" || echo "$PWD/${1#./}"; }
usage() {
  die "Usage: $0 -n [target sequences] -w [work directory] -r [check for reqs]."
}
exit_abnormal() {
  usage
  exit 1
}

OPTIONS=$(getopt -o :n:w:hr -l number:,workdir:,help,reqs -- "$@")


if [ $? -ne 0 ]; then
  usage
  exit 1
fi

eval set -- $OPTIONS

while true; do
  case "$1" in

    -w | --workdir)
    WORKDIR="$(realpath $2)"
    shift
    ;;
    -n | --number)
    MAX_TARGET_SEQS="$2"
    shift
    ;;
    -r | --reqs)

    probe git || die
    probe wget || die
    probe dos2unix || die
    probe blastn || die
    probe timeout || die
    probe jq || die

    exit 0;
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
    echo $OPTIONS
    exit 1
    ;;
  esac
  shift
done

probe git
probe wget
probe dos2unix
probe blastn
probe timeout
probe jq

INPUT="fasta";

if [ ! -z "$@" ]; then
  STATUS="Warning :"
  yell "\"$@\" is extraneous. Ignoring."
fi

if [ -z $WORKDIR ]; then

  WORKDIR=$(realpath .)

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

  yell "Defaulting maximum target sequences to ${MAX_TARGET_SEQS}, was ${invalid_value}."

else
  yell "MAX_TARGET_SEQS is ${MAX_TARGET_SEQS}."

fi

#Make sure all commands work
STATUS="Probe:"
BLASTPATH="$(realpath "$0")"
BLASTDIR="$(dirname "$BLASTPATH")"
ARGS="$@"

body() {
  IFS= read -r header
  printf '%s\n' "$header"
  "$@"
}

cd "$WORKDIR";
mkdir quickblast_tmp 2>/dev/null
cd quickblast_tmp || die
WORKDIR=$(realpath .)

awk 'BEGIN{RS=">";count=1} $0~/[A-z]/{print ">"$0 > "seq_"count".fa";count++}' ../*.fa* 2>/dev/null

  STEP="BlastN:"
  PROGRESS=0;

  COUNT=$(ls -1 *.fa* 2>/dev/null | wc -l | tr -d ' ')

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
          timeout --foreground 5m blastn -db nt -query $x -remote -max_target_seqs=${MAX_TARGET_SEQS} -out $OUTFILE -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" 2>&1 |

          while read line; do
            STATUS="NCBI returned:"
            yell $line

            if [[ $line == *"failed"* ]] || [[ $line == *"Failed"* ]]; then
              ps -ef | grep $x | grep $USER | grep timeout | grep -v grep | awk '{print $2}' | xargs kill && break

            fi


          done

          STATUS="($PROGRESS/$COUNT) Result:"
          newOUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')
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
          timeout --foreground 5m blastn -db nt -query $x -remote -max_target_seqs=${MAX_TARGET_SEQS} -out $OUTFILE -outfmt "6 qseqid stitle sacc sseqid pident qlen length evalue bitscore" 2>&1 |

          while read line; do
            STATUS="($PROGRESS/$COUNT) NCBI returned:"
            yell $line

            if [[ $line == *"failed"* ]] || [[ $line == *"Failed"* ]]; then
              ps -ef | grep $x | grep $USER | grep timeout | grep -v grep | awk '{print $2}' | xargs kill && break

            fi

          done

          STATUS="($PROGRESS/$COUNT) Result:"
          newOUTSIZE=$(cat $OUTFILE | wc -l | tr -d ' ')
          yell "${OUTFILE}, size: ${newOUTSIZE} lines."
          yell "Pausing..."
          sleep 60

          STATUS="($PROGRESS/$COUNT) Blasting:"


          if [ "$newOUTSIZE" = "0" ]; then
            STATUS="($PROGRESS/$COUNT) Warning:"
            yell "Unable to initialize ${OUTFILE}"
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
  mv blast_out.csv ../
  yell "Blast results located at \"blast_out.csv\"."

fi
