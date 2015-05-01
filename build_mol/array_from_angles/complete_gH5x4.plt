set terminal png enhanced;set output "complete_gH5x4.png";
set ylabel '';set format y '';set xtics nomirror;set ytics nomirror; set border 3
set style line 11 lc rgb '#808080' lt 1;set border 3 back ls 11;
plot 'complete_gH5x4.pdb.dat' u 1:2 thru log(y) t 'FoXS' w lines lw 2.5 lc rgb '#e26261'
reset
