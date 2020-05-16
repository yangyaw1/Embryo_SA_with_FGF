1: The name for the gene file is "XXX_gene.timestamp.pqr" in "signal" folder where XXX is the file name you can modify in main file. The timestamp is the simulation time (16C: 39000, 32C: 75000, 64C: 95000, 128C beginning: 120000, 128C mid: 160000, end of simulation: 260000)

2: The file contains several rows and columns. Each row is one cell. The column shows gene expression information: 
col 1: cell ID
col 2: cell type (On or after 16C: 0: EPI, 1: DP, 2: TE, 4: PE; 
                  1C-8C, cell type is 1 representing 'undifferentiated' cell)
col 3: Cdx2 expression
col 4: Oct4 expression
col 5: Nanog expression
col 6: Gata6 expression
col 7: Fgfr2 expression
col 8: Erk expression
col 9: Fgf4 expression inside the cell
col 10: Fgf4 expression received from the neighborhood 

