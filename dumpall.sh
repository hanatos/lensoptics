#!/bin/bash
degree=11
maxcoeff=40
./genpoly $degree sorted.poly
suffix=""
lenses=$1
if [ "$1" == "" ]
then
  lenses=lenses/*.fx
fi
for i in $lenses
do
  lens1=${i#*/}
  lens=${lens1%.*}
  echo '[' $lens ']'
  ./fit $i $degree $degree $degree $maxcoeff
  ./gencode $i
  ./view $i -o
  ./fresnel $i
  gnuplot -e "lens=\"$i\"" fresnel.plt
  mv screenshot.pdf lenses/${lens}.pdf
  mv fresnel-profile.pdf lenses/${lens}-zzz-fresnel.pdf
  echo ""
  mkdir -p render/${lens}${suffix}
  mv *.h render/${lens}${suffix}/
done
pdftk lenses/*pdf cat output lenses.pdf

