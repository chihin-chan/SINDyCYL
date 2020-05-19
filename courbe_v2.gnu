set terminal png

list = system('ls vort0*')
i = 1

do for [file in list] {
    #set output sprintf('v0/%s.png', file)
    set output 'res_'.file.'.png'
    #set output sprintf('_/%s.png', file)
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
    set size 1,1
    splot file u 1:2:3
    reset
    i = i + 1
}
