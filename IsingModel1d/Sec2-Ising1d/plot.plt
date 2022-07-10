set terminal gif animate delay 2
set output 'movie.gif'
set xrange [0:L]
set yrange [-1:1]
set size ratio 0.2
set xlab "lattice site"
set ylab "magnetisation"
set ytics 1

set palette defined ( -1 '#000fff', 1 '#ee0000')
unset colorbox

do for [i=0:10000:10] {
    plot 'out.'.i.'.dat' u 1:(-$2/2):(0):2:2 w vec palette lw 0.2 not
}
