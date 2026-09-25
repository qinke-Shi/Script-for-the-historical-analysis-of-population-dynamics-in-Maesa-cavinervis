// 5-group model M21: Bottleneck in G5 and G6 + Gene flow between G2 and G3 + Gene flow between G5 and G6
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
3
//Migration matrix 0 (gene flow: G2<->G3 and G5<->G6)
0 MIG_G6_to_G5 0 0 0
MIG_G5_to_G6 0 0 0 0
0 0 0 MIG_G2_to_G3 0
0 0 MIG_G3_to_G2 0 0
0 0 0 0 0
//Migration matrix 1 (gene flow: G2<->G3 only; used after G6 merges into G5 at TIME1)
0 0 0 0 0
0 0 0 0 0
0 0 0 MIG_G2_to_G3 0
0 0 MIG_G3_to_G2 0 0
0 0 0 0 0
//Migration matrix 2 (no gene flow; used after G2 merges into G1 at TIME3)
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
6 historical event
TBOT_E 0 0 0 RES_G5 0 0
TBOT_E 1 1 0 RES_G6 0 0
TIME1 1 0 1 1 0 1
TIME2 0 2 1 1 0 1
TIME3 3 4 1 1 0 2
TIME4 2 4 1 1 0 2
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block
1
//per Block:data type, number of loci, recomb and mut rates
FREQ 1 0 2.43e-8 OUTEXP
