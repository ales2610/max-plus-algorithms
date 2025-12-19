function [delta_min, delta_max, low_95, upp_95] = lea(A, c, j, N_max, M, delta, pr)
% LEA Monte Carlo estimation of transient growth bounds
% in max-plus linear systems.
%
% This function estimates lower and upper bounds on the deviation
% between two max-plus trajectories:
% - one starting at time 0
% - one starting at time c
%
% The method repeatedly perturbs the system matrix A, propagates
% max-plus products, and measures the difference in the j-th column.
%
% INPUT:
% A : n x n max-plus matrix (entries in R ∪ {-Inf})
% c : time shift (delay) parameter
% j : reference column index
% N_max : maximum time horizon
% M : number of Monte Carlo samples
% delta : stopping tolerance on (delta_max - delta_min)
% pr : plotting flag (true/false)
%
% OUTPUT:
% delta_min : estimated lower bound on average deviation
% delta_max : estimated upper bound on average deviation
% low_95 : 2.5%% percentile of lower deviation estimates
% upp_95 : 97.5%% percentile of lower deviation estimates
%
% The algorithm is related to Lyapunov-type and coupling-time
% analysis for max-plus stochastic systems.

    n = size(A,1);

    % Accumulators for averaged bounds
    delta_min = 0;
    delta_max = 0;

    % Average evolution of bounds over time
    avg_d = zeros(2,N_max);
    % Store all Monte Carlo samples
    d_all = zeros(2,M);

    % --- Monte Carlo loop --
    for m=1:M
        % Bounds for the m-th sample
        delta_m_min = -inf;
        delta_m_max = inf;

        % Initialize max-plus identity matrices
        A_k0 = identity(n);             % product starting at time 0
        A_kc = identity(n);             % product starting at time c

        % --- Time evolution ---
        for k=0:N_max-1
          
            % Apply random perturbation to A
            A_k = exponential(A);          

            % Accumulate products
            A_k0 = otimes(A_k, A_k0);
            
            % Accumulate delayed product only after time c
            if k >= c
                A_kc = otimes(A_k, A_kc);
            end

            % Difference between the two trajectories
            Delta = A_k0 - A_kc;

            % Lower bound from column j
            delta_m_min = min(Delta(:,j));
            if isnan(delta_m_min), delta_m_min = -inf; end
     
            % Upper bound from column j
            delta_m_max = max(Delta(:,j));
            if isnan(delta_m_max), delta_m_max = inf; end

            % Accumulate averages
            avg_d(1, k+1) = avg_d(1, k+1) + delta_m_min;
            avg_d(2, k+1) = avg_d(2, k+1) + delta_m_max;

            % Early stopping if bounds are tight enough
            if ~pr && delta_m_max - delta_m_min <= delta, break; end
        end

        % Accumulate Monte Carlo averages
        delta_min = delta_min + delta_m_min;
        delta_max = delta_max + delta_m_max;

        % Store individual samples
        d_all(1,m) = delta_m_min;
        d_all(2,m) = delta_m_max;
    end

    % Normalize by time shift and number of samples
    delta_min = delta_min / (c * M);
    delta_max = delta_max / (c * M);

    avg_d = avg_d ./ (c * M);

    d_all = d_all ./ c;

    % 95% confidence interval (lower bound distribution)
    low_95 = prctile(d_all(1,:), 2.5);
    upp_95 = prctile(d_all(1,:), 97.5);

    % First index where bounds become finite
    idx = find(isinf(avg_d(1, :)), 1, 'first');

    % Optional plotting
    if pr
        x = idx:1:N_max;

        figure
        plot(x, avg_d(1, idx:end), x, avg_d(2, idx:end))

        figure
        histogram(d_all(1,1:end));
    end
end


% -------------------------------------------------------------------------
function E = identity(n)
% IDENTITY Max-plus identity matrix
% E(i,i) = 0, E(i,j) = -Inf for i ~= j


E = -inf*ones(n);
for i=1:n
E(i,i) = 0;
end
end


% -------------------------------------------------------------------------
function C = otimes(A, B)
% OTIMES Max-plus matrix multiplication
% C = A ⊗ B, where C(i,j) = max_k (A(i,k) + B(k,j))


n = size(A,1);
C = -inf(n);


for k = 1:n
C = max(C, A(:,k) + B(k,:));
end
end


% -------------------------------------------------------------------------
function F = exponential(A)
% EXPONENTIAL Random exponential perturbation of matrix entries
%
% Each finite entry A(i,j) is randomly perturbed according to
% an exponential distribution, modelling stochastic processing times.


p = 0.7; % probability of keeping an entry unchanged
mu = 0.05; % exponential rate parameter


F = A;
R = rand(size(A));
U = rand(size(A));


% Apply perturbation only to finite entries
mask = (R > p) & ~isinf(A);


% Exponential noise
F(mask) = A(mask) + (-(mu .* A(mask)) .* log(U(mask)));
end