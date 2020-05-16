1: The name for the space file is "XXX_space.timestamp.pqr" in "PQR" folder where XXX is the file name you can modify in main file. The timestamp is the simulation time (16C: 39000, 32C: 75000, 64C: 95000, 128C beginning: 120000, 128C mid: 160000, end of simulation: 260000)

2: The file contains several rows and columns. Each row is one element. The column shows cell and space information: 
col 2: the element id
col 5: the mixture of cell id and cell type, [the num of col]%200 is the cell id,
	[the num of col]/1000 is the cell type (1 is Nanog+ cell, 2 is DP cell,
						3 is TE cell, 5 is Gata+ cell)
col 6, 7, 8: the x, y, z coordinate is 3D space
col 9, 10: related to element type but not used in the code (needed when you use VMD to 		draw the embryo)

3: You can use VMD to draw the embryo with the help of "drow_pqr.sh". 

