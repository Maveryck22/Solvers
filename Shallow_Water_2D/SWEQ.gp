reset
set terminal pngcairo size 900,524 enhanced font 'Verdana,10'

# initializing values for the loop and start the loop

t = 2
end_time = 3
system('mkdir -p animation')
load 'SWEQ2D.plt'
