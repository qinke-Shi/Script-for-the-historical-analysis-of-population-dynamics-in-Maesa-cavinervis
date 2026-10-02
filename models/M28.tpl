// 5-group model M28: Bottleneck in G5 and G6, no gene flow
5 samples to simulate :
//Population effective sizes (number of genes)
G5N
G6N
G3N
G2N
G1N
//Samples sizes and samples age
18
50
22
34
28
//Growth rates
0
0
0
0
0
//Number of migration matrices
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
8 historical event
TBOT_G5 0 0 0 botr_G5 0 0
TENDBOT_G5 0 0 0 recr_G5 0 0
TBOT_G6 1 1 0 botr_G6 0 0
TENDBOT_G6 1 1 0 recr_G6 0 0
TIME1 1 0 1 1 0 0
TIME2 0 2 1 1 0 0
TIME3 3 4 1 1 0 0
TIME4 2 4 1 1 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block
1
//per Block:data type, number of loci, recomb and mut rates
FREQ 1 0 2.43e-8 OUTEXP
