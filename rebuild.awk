#!/usr/bin/awk -E
BEGIN {
  FS="\t";
  OFS=","
} {
  rebuilt=0;
  for(i=1; i<=NF; ++i) {
    if ($i ~ /,/ && $i !~ /^".*"$/) {
      $i = "\"" $i "\"";
      rebuilt=1
    }
  }
  if (!rebuilt) {
    $1=$1
  }
  print
}
