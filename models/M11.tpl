// 5-group model M11: Bottleneck in G5 and G6 + Gene flow between G5(0) and G3(2)
5 samples to simulate :
//Population effective sizes (number of genes)
G5N
G6N
G3N
G2N
G1N
//Samples sizes and samples age
22
50
22
34
30
//Growth rates
0
0
0
0
0
//Number of migration matrices
2
//Migration matrix 0 (gene flow between G5(0) and G3(2))
0 0 MIG_G3_to_G5 0 0
0 0 0 0 0
MIG_G5_to_G3 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//Migration matrix 1 (no gene flow)
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
6 historical event
TBOT_E 0 0 0 RES_G5 0 0
TBOT_E 1 1 0 RES_G6 0 0
TIME1 1 0 1 1 0 0
TIME2 0 2 1 1 0 1
TIME3 3 4 1 1 0 1
TIME4 2 4 1 1 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block
1
//per Block:data type, number of loci, recomb and mut rates
FREQ 1 0 2.43e-8 OUTEXP
