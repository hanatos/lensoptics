set terminal png
set output 'fresnel-rt.png'
set size ratio 1
#set palette 1 1
#set xrange [-18:18]
#set yrange [-12:12]
set key bottom right
plot 'fresnel.dat' u 1:2:3 w p pt 5 lw 3 ps 0.01 lc palette title 'transmittance rt'
set output 'fresnel-poly.png'
plot 'fresnel.dat' u 4:5:6 w p pt 5 lw 3 ps 0.01 lc palette title 'transmittance poly'
set terminal pdf
set title 'radial transmittance profile'
set output 'fresnel-profile.pdf'
plot 'fresnel-profile.dat' u 0:1 w l title 'ray traced',\
     'fresnel-profile.dat' u 0:2 w l title 'polynomial'
unset output
#!pdfcrop 'fresnel.pdf' 'fresnel.pdf'
