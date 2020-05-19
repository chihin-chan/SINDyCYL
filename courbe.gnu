set terminal png
set output 'v01.png'
set xr[0:20] #0.025
set yr[0:12]  #0.15
set cbrange [-0.25:0.25]
unset colorbox
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
#set grid lw 2 lt 6
set xtics out
set ytics out
set pm3d map
#set pm3d interpolate 0,0
set size 1.10,1
splot 'vort0001' u 1:2:3

reset
set terminal png
set output 'v02.png'
set xr[0:20] #0.025
set yr[0:12]  #0.15
set cbrange [-0.25:0.25]
unset colorbox
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
#set grid lw 2 lt 6
set xtics out
set ytics out
set pm3d map
#set pm3d interpolate 0,0
set size 1.1,1
splot 'vort0002' u 1:2:3
reset

set terminal png
set output 'v03.png'
set xr[0:20] #0.025
set yr[0:12]  #0.15
set cbrange [-0.25:0.25]
unset colorbox
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
#set grid lw 2 lt 6
set xtics out
set ytics out
set pm3d map
#set pm3d interpolate 0,0
set size 1.1,1
splot 'vort0003' u 1:2:3

reset
set terminal png
set output 'v04.png'
set xr[0:20] #0.025
set yr[0:12]  #0.15
set cbrange [-0.25:0.25]
unset colorbox
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
#set grid lw 2 lt 6
set xtics out
set ytics out
set pm3d map
#set pm3d interpolate 0,0
set size 1.1,1
splot 'vort0004' u 1:2:3