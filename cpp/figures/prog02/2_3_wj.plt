#set term qt persist
set term tikz color solid size 5in,3in
set output "2_3_wj.tikz"

set title "2.3 - Weighted Jacobi"

set key

set xlabel "Iterations"
set ylabel "Error - $\\|\\mathbf{e}\\|_\\infty$" noenhanced

plot "2_3_wj_k1_error.txt" using 1:2 with lines\
         linetype rgb "black" title "$k = 1$",\
     "2_3_wj_k3_error.txt" using 1:2 with lines\
         linetype rgb "red" title "$k = 3$",\
     "2_3_wj_k6_error.txt" using 1:2 with lines\
         linetype rgb "blue" title "$k = 6$"

