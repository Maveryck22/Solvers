# bessel loop
t = t + 100

outfile = sprintf('animation/bessel%02.0f.png',t)
set output outfile

# --------------- archivo comun para plotear ------------------

set dgrid3d 20,20
set pm3d
set hidden3d offset 0 # Se agregan las lineas de malla
set ticslevel 0.0
set isosample 15,40

set xlabel 'x'
set ylabel 'y'
set zlabel 'V'

show dgrid3d
set view 60,135*0.3*t/20000,1 # Giros eje horizontal, eje z, zoom

set xrange[0:1]
set yrange[0:1]
set zrange[0.4:1]

set palette defined (-1 "cyan", -0.5 "blue", 0 "green",0.4 "yellow",1 "red")

splot   'VelocidadU.dat' u 1:2:t w pm3d ls 1 title 'Burgers Equation'


if(t<end_time) reread;
