set terminal png
set output 'fresnel.png'
set size ratio 1
#set palette 1 1
#set xrange [-18:18]
#set yrange [-12:12]
set key bottom right
plot 'fresnel.dat' w p pt 5 lw 3 ps 0.01 lc palette title 'transmittance'
unset output
#!pdfcrop 'fresnel.pdf' 'fresnel.pdf'
