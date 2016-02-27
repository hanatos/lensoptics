#!/bin/bash
degree=5
suffix=""
for i in lenses/*.fx
do
  lens1=${i#*/}
  lens=${lens1%.*}
  echo '[' $lens ']'
  ./fit $i $degree
  ./fit $i -2
  ./simplify $i
  ./fit $i -2
  ./gencode $i
  ./view $i -o
  mv screenshot.pdf lenses/${lens}.pdf
  echo ""
  mkdir -p render/${lens}${suffix}
  mv *.h render/${lens}${suffix}/
done
pdftk lenses/*pdf cat output lenses.pdf

