function [eta, v, j] = howard(A)
% HOWARD  Implementation of Howard's policy iteration algorithm
% for computing a generalised max-plus eigenmode.
%
% This function computes:
%   - eta : the max-plus eigenvalue (cycle mean) associated with each state
%   - v   : the eigenvector
%   - j   : an index of a critical node
%
% INPUT:
%   A : 3D array of size (n x n x m)
%       A(i,j,l) is the weight of the transition j -> i with delay l.
%       Non-existing transitions are -Inf.
%
% OUTPUT:
%   eta : n x 1 vector of eigenvalues
%   v   : n x 1 vector of eigenvector values
%   j   : index of a critical node
%
% The algorithm alternates between:
%   1) Policy evaluation (value determination)
%   2) Policy improvement (Howard iteration)


    %start clock
    tic

    % Problem dimensions
    n = size(A,1);      % number of states
    m = size(A,3);      % order of system

    % Initial value function
    v = zeros(n,1);

    % A_pi stores the transition matrix induced by the current policy pi
    A_pi = -inf * ones(n, n, m);

    % Policy pi(i,:) = [next_state, duration_index]
    pi = zeros(n, 2);

    % Initialize policy: choose the first finite transition for each state
    for i = 1:n
        [j_idx, l_idx] = ind2sub([n, m], find(isfinite(A(i,:,:))));
        pi(i, 1) = j_idx(1);
        pi(i, 2) = l_idx(1);
        A_pi(i, pi(i, 1), pi(i,2)) = A(i, pi(i, 1), pi(i, 2));
    end
    
    % Perform Howard policy iteration
    [eta, v, critical_j] = improvement(A, pi, A_pi, v);
    j = critical_j;

    %stop clock
    elapsed = toc;
    fprintf('Time of Howard: %.6f seconds\n', elapsed);
end

% -------------------------------------------------------------------------
function xi = cycle(A, s)
% CYCLE  Detects a cycle reachable from state s in communication graph of matirx A

    n = size(A,1);
    m = size(A,3);
    visited = false(1,n);
    
    % xi stores the visited states and actions
    xi = zeros(n+1,2);
    xi(1,1) = s;
    visited(s) = true;
    
    k = 1;
    while true
        % Extract outgoing transition from current state
        A_i = squeeze(A(xi(k),:,:))';
        [l, j] = ind2sub([m n], find(isfinite(A_i)));
        
        % Store chosen action and next state
        xi(k, 2) = l;
        xi(k+1, 1) = j;
    
        % If next state was already visited, a cycle is found
        if visited(xi(k+1,1))
            cycle_start = find(xi(:,1) == xi(k+1,1), 1);
            xi = xi(cycle_start:k+1, :);
            break
        end
    
        visited(xi(k+1, 1)) = true;
        k = k + 1;
    end
end

% -------------------------------------------------------------------------
function [eta, v, visited] = visit_all(A_pi, xi, eta, v, j, visited)
% VISIT_ALL  Propagates eigenvalue and eigenvector values
% from a critical node to all reachable nodes.

    n = size(A_pi,1);
    m = size(A_pi,3);

    % Transitions going into node j
    A_j = squeeze(A_pi(:,j,:));
    [idx_i, idx_l] = ind2sub([n m], find(isfinite(A_j)));
        
    for k = 1:numel(idx_i)
        i = idx_i(k);
    
        if ~visited(i)
            visited(i) = 1;

            % Same eigenvalue as the successor
            eta(i) = eta(j);

            % Update for the eigenvector
            v(i) = A_pi(i,j,idx_l(k)) - ((idx_l(k) - 1) * eta(j)) + v(j);
        
            % Recursive propagation
            [eta, v, visited] = visit_all(A_pi, xi, eta, v, i, visited);
        end
    end
end

% -------------------------------------------------------------------------
function eta_bar = avg_weight(xi, A_pi)
% AVG_WEIGHT  Computes the average weight (cycle mean)
% of a critical cycle.

    eta_bar = 0;
    len = 0;

    for i = 1:length(xi(:,1))-1
        eta_bar = eta_bar + A_pi(xi(i,1), xi(i+1,1), xi(i,2));

        %lenght is defined as sum of all delays
        len = len + (xi(i,2) - 1);
    end
    
    eta_bar = eta_bar / len;
end

% -------------------------------------------------------------------------
function [eta, v, j] = value_det(A_pi, v_old)
% VALUE_DET  Policy evaluation step.
%
% Computes eigenvalues and eigenvectors for the current policy
% by identifying critical cycles and propagating values.

    n = size(A_pi,1);

    eta = -inf * ones(n,1);
    v = -inf * ones(n,1);

    j = 0;  % representative critical node

    s = 1;  % starting state
    visited = zeros(1, n);

    while ~isempty(s)
        % Find a cycle reachable from s
        xi = cycle(A_pi, s);

        % Compute cycle mean
        eta_bar = avg_weight(xi, A_pi);

        % Store a critical node
        if j == 0, j = xi(1); end
        
        % Initialize one node on the cycle
        eta(xi(1,1)) = eta_bar;
        v(xi(1,1)) = v_old(xi(1,1));
        
        visited(xi(1,1)) = 1;

        % Propagate values to all reachable nodes
        [eta, v, visited] = visit_all(A_pi, xi(1:end-1,:), eta, v, xi(1,1), visited);

        % Continue with an unvisited node
        s = find(visited == 0, 1);
    end
end

% -------------------------------------------------------------------------
function [eta, v, j] = improvement(A, pi, A_pi, v)
% IMPROVEMENT  Policy improvement loop of Howard's algorithm.
%
% Alternates between:
%   - improving the eigenvalue (eta)
%   - improving the eigenvector (v)
% until no further improvement is possible.

    tol = 1e-10;

    n = size(A_pi,1);
    m = size(A_pi,3);
    cont_improvement = true;

    while cont_improvement

        cont_improvement1 = true;

        % --- Eigenvalue (eta) improvement ---
        while cont_improvement1
            [eta, v, j] = value_det(A_pi, v);

            I_1 = 0;
            for i=1:n
                A_i = squeeze(A(i,:,:))';
                eta_j = repmat(eta', m, 1);

                eta_j(~isfinite(A_i)) = -inf;

                [max_eta, idx] = max(eta_j(:));
                [max_l, max_j] = ind2sub([m n], idx);
    
                % Find nodes that do not match the condidtion
                if eta(i) < max_eta - tol
                    I_1 = I_1 + 1; 

                    % Update policy
                    A_pi(i,pi(i,1),pi(i,2)) = -inf;
                    pi(i,1) = max_j;
                    pi(i,2) = max_l;
                    A_pi(i,max_j,max_l) = A(i,max_j,max_l);
                end
            end
    
            if I_1 == 0
                % P_bar: transitions preserving maximal eta
                P_bar = (eta(:) == eta(:).') .* isfinite(A);
                cont_improvement1 = false;
            end    
        end

        % --- Eigenvector (v) improvement ---
        I_2 = 0;
        for i=1:n
            A_i = squeeze(A(i,:,:))';
            P_bar_i = squeeze(P_bar(i,:,:))';
            v_j = repmat(v', m, 1);
            eta_j = repmat(eta', m, 1);

            val_i = A_i + v_j - (0:(m-1))' .* eta_j;
            val_i(P_bar_i == 0) = -inf;
            
            [max_v, idx] = max(val_i(:));
            [max_l, max_j] = ind2sub([m n], idx);
            
            % Find nodes that do not match the condidtion
            if v(i) < max_v - tol
                I_2 = I_2 + 1; 

                % Update policy
                A_pi(i,pi(i,1),pi(i,2)) = -inf;
                pi(i,1) = max_j;
                pi(i,2) = max_l;
                A_pi(i,max_j,max_l) = A(i,max_j,max_l);
            end
        end
    
        if I_2 == 0
            cont_improvement = false;
        end     
    end
end
