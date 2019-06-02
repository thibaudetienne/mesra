#!/bin/bash

nstates=`grep -i 'nstates' $1`
nstates_nws="$(echo -e "${nstates}" | tr -d '[:space:]')"

nsvar=`echo -e $nstates_nws | tr [:upper:] [:lower:] | awk '{split($0,a,"nstates=") ; print a[2]}'`

for (( i=0; i<${#nsvar}; i++ )); do
  echo "${nsvar:$i:1}" >> ns_letter_by_letter
done

for f in $(cat ns_letter_by_letter)
do
 if [[ $f == "," ]];then
  exit
 elif [[ $f == " " ]];then
  exit
 elif [[ $f == ")" ]];then
  exit
 else echo $f
 fi
done

rm ns_letter_by_letter

