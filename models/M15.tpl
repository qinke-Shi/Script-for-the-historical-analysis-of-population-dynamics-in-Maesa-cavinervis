// 5-group model M15: No bottleneck + Gene flow between G2 and G3 (like M03) + Gene flow between G1 and G3
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
//Migration matrix 0 (gene flow: G2<->G3 and G1<->G3)
0 0 0 0 0
0 0 0 0 0
0 0 0 MIG_G2_to_G3 MIG_G1_to_G3
0 0 MIG_G3_to_G2 0 0
0 0 MIG_G3_to_G1 0 0
//Migration matrix 1 (no gene flow)
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
4 historical event
TIME1 1 0 1 1 0 0
TIME2 0 2 1 1 0 0
TIME3 3 4 1 1 0 1
TIME4 2 4 1 1 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block
1
//per Block:data type, number of loci, recomb and mut rates
FREQ 1 0 2.43e-8 OUTEXP
