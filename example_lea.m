% Initialize matrix
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

% Parameters of this example
c = 4;
j = 3;
delta = 10^(-8);
N_max = 250;
M = 1000;

% Use LEA algorithm to get lower and upper bound for Lyapunov exponent. For this example both should be around 53.66345.
[delta_min, delta_max, low_95, upp_95] = lea(A, c, j, N_max, M, delta, true)

