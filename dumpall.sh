#!/bin/bash
# suffix="-taylor"
suffix="-fit"
for i in lenses/*.fx
do
  lens1=${i#*/}
  lens=${lens1%.*}
  echo '[' $lens ']'
# ./fit $i
#  ./dump-code $i >& /dev/null
  ./view $i -o
  mv screenshot.pdf lenses/${lens}.pdf
  echo ""
# mkdir -p ~/vcs/corona-6/camera/${lens}${suffix}
#  mv *.h ~/vcs/corona-6/camera/${lens}${suffix}/
done
pdftk lenses/*pdf cat output lenses.pdf

