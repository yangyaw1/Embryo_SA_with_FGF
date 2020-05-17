# Embryo development simulation with 3D model from zygote to late blastocyst
A study of early mammalian embryo development using subcellular element method for spatial dynamics and gene network ODE models for gene expression dynamics.
The two main aspects we focus on are 1) selective adhesion and 2) FGF signaling.

## How to execute the code
1: The code is executed using <code>C</code> with <code>OpenCL</code>.

2: To execute the code, please copy the following into the commend line:
```
gcc -o exe main.c -lOpenCL -std=c99 -lm
./exe
```

3: You need to choose the platform to execute the code in line 474 of <code>main.c</code>. In the same line, use <code>CL_DEVICE_TYPE_GPU</code> for GPU device or <code>CL_DEVICE_TYPE_CPU</code> for CPU device.
(It is show as platform_id[0] in the code, if you have error when executing,
you can try platform_id[1], platform_id[2], ...)

4: You can set the parameter related to the results of the paper at line 20-line 44 of <code>main.c</code>

5: Once you executed the code, you will get space files (in "PQR" folder) and gene files (in "signal" folder).

## How to use the space file

1: The name for the space file is <code>XXX_space.timestamp.pqr</code> in <code>PQR</code> folder where XXX is the file name you can modify in main file. The timestamp is the simulation time (16C: 39000, 32C: 75000, 64C: 95000, 128C beginning: 120000, 128C mid: 160000, end of simulation: 260000)

2: The file contains several rows and columns. Each row is one element. The column shows cell and space information:  

col 2: the element id  

col 5: the mixture of cell id and cell type, [the num of col]%200 is the cell id, [the num of col]/1000 is the cell type (1 is Nanog+ cell, 2 is DP cell, 3 is TE cell, 5 is Gata+ cell)  

col 6, 7, 8: the x, y, z coordinate is 3D space  

col 9, 10: related to element type but not used in the code (needed when you use VMD to draw the embryo)

3: You can use <code>VMD</code> to draw the embryo. First load the pqr file as new molecule into <code>VMD</code>, then execute the commands of <code>draw_pqr.sh</code> in Tk console. 

## How to use the gene file

1: The name for the gene file is <code>XXX_gene.timestamp.txt</code> in <code>signal</code> folder where XXX is the file name you can modify in main file. The timestamp is the simulation time (16C: 39000, 32C: 75000, 64C: 95000, 128C beginning: 120000, 128C mid: 160000, end of simulation: 260000)

2: The file contains several rows and columns. Each row is one cell. The column shows gene expression information:  

col 1: cell ID  

col 2: cell type (On or after 16C: 0: EPI, 1: DP, 2: TE, 4: PE;  1C-8C, cell type is 1 representing 'undifferentiated' cell)  

col 3: Cdx2 expression  

col 4: Oct4 expression  

col 5: Nanog expression  

col 6: Gata6 expression  

col 7: Fgfr2 expression  

col 8: Erk expression  

col 9: Fgf4 expression inside the cell  

col 10: Fgf4 expression received from the neighborhood  
