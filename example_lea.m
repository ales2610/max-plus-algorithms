%initialize matrix
A = [-inf -inf -inf -inf -inf  74  -inf  78  -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf  210 -inf  192 -inf
     -inf  48  -inf  61  -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf 
     -inf -inf -inf -inf -inf -inf -inf -inf  160 -inf  148 -inf -inf -inf -inf -inf 
       0  -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf
     -inf   0  -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf   0  -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf   0  -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf   0  -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf   0  -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf -inf   0  -inf -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf -inf -inf   0  -inf -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf -inf -inf -inf   0  -inf -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf -inf -inf -inf -inf   0  -inf -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf   0  -inf -inf -inf -inf -inf
     -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf   0  -inf -inf -inf -inf];

%parameters of this example
c = 4;
j = 3;
delta = 10^(-8);
N_max = 250;
M = 1000;

[delta_min, delta_max, low_95, upp_95] = lea(A, c, j, N_max, M, delta, false)
