Use "#drawe.sh" file to draw Nanog+,Gata6+,DP
Use "#drawe_fgf.sh" file to draw Fgf4+/-,Fgfr2+/-

Useful information are in columns 5, 6, 7, 8

The cell type determined in column 5:
id > 5000: Gata6+
3000 < id <= 5000: PE
2000 < id <= 3000: DP
else: Nanog+

The fgf4/fgfr2 info determined in column 5:
for id which is not PE:
val = id%1000
if val >= 600: fgf4+/fgfr2+
if 400 <= val < 600: fgf4-/fgfr2+
if 200 <= val < 400: fgf4+/fgfr2-


For figure 2d and 3b
h1, h2, h7 in fig2d only contain 5 simulations since they must be fail. The others contain 10 simulations.
 