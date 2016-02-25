#!/bin/bash
suffix=""
for i in lenses/*.fx
do
  lens1=${i#*/}
  lens=${lens1%.*}
  echo '[' $lens ']'
  ./fit $i
  ./gencode $i
  ./view $i -o
  mv screenshot.pdf lenses/${lens}.pdf
  echo ""
  mkdir -p ../corona-13/camera/${lens}${suffix}
  mv *.h ../corona-13/camera/${lens}${suffix}/
done
pdftk lenses/*pdf cat output lenses.pdf

