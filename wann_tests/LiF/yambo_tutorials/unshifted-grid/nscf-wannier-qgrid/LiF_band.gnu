set style data dots
set nokey
set xrange [0: 4.73667]
set yrange [-20.82450 : 23.71879]
set arrow from  1.08988, -20.82450 to  1.08988,  23.71879 nohead
set arrow from  2.42470, -20.82450 to  2.42470,  23.71879 nohead
set arrow from  3.96601, -20.82450 to  3.96601,  23.71879 nohead
set xtics ("W"  0.00000,"L"  1.08988,"Γ"  2.42470,"X"  3.96601,"W"  4.73667)
 plot "LiF_band.dat"
