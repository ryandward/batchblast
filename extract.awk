#!/usr/bin/awk -E

BEGIN{
  OFS=",";
  FPAT = "([^,]+)|(\"[^\"]+\")"
}
{
  tax = ""

  while (tax==""){
    curl[NR]="curl -s http://taxonomy.jgi-psf.org/tax/accession/"$3" | jq .[][.[].level].name"
    curl[NR] | getline tax;
    close(curl[NR]);

  }
  $13=tax
  $10=$9
  $9=$8
  $8=$7
  $7=$6
  $6=$5
  $11="=HYPERLINK(\"https://www.ncbi.nlm.nih.gov/nuccore/"$3"?report=fasta\")";
  $12=$2
  $2=$3
  gsub(/"/, "",tax)
  gsub(/'/, "",tax)
  split(tax,parts," ")

  $3=parts[1]

  $5=""

  if(parts[2]=="cf."){
    $4="\""parts[2]" "parts[3]"\""
    gsub(parts[1],"",tax)
    gsub(parts[2],"",tax)
    gsub(parts[3],"",tax)
    gsub(/^[[:space:]]+/,"",tax)
    $5=tax

  }

  else{
    $4=parts[2]
    gsub(parts[1],"",tax)
    gsub(parts[2],"",tax)
    gsub(/^[[:space:]]+/,"",tax)
    $5=tax

  }

  if(index($5, " ")){
    $5="\""$5"\""
  }

  for(i in parts){
    delete parts[i]
  }
  print $0;
  tax=0;
  next;
}
