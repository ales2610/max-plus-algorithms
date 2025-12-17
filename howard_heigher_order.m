function [eta, v] = howard_heigher_order(A)

    n = size(A,1);
    m = size(A,3);
    v = zeros(n,1);

    A_pi = -inf * ones(n, n, m);
    pi = zeros(n, 2);

    for i = 1:n
        [j, l] = ind2sub([n, m], find(isfinite(A(i,:,:))));
        pi(i, 1) = j(1);
        pi(i, 2) = l(1);
        A_pi(i, pi(i, 1), pi(i,2)) = A(i, pi(i, 1), pi(i, 2));
    end
    
    [eta, v, critical_j] = improvement(A, pi, A_pi, v);
end

function xi = cycle(A, s)

    n = size(A,1);
    m = size(A,3);
    visited = false(1,n);
    
    xi = zeros(n+1,2);
    xi(1,1) = s;
    visited(s) = true;
    
    k = 1;
    while true
        A_i = squeeze(A(xi(k),:,:))';
        [l, j] = ind2sub([m n], find(isfinite(A_i)));
        
        xi(k, 2) = l;
        xi(k+1, 1) = j;
    
        if visited(xi(k+1))
            cycle_start = find(xi(:,1) == xi(k+1,1), 1);
            xi = xi(cycle_start:k+1, :);
            break
        end
    
        visited(xi(k+1, 1)) = true;
        k = k + 1;
    end
end

function [eta, v, visited] = visit_all(A_pi, xi, eta, v, j, visited)

    n = size(A_pi,1);
    m = size(A_pi,3);

    A_j = squeeze(A_pi(:,j,:));
    [idx_i, idx_l] = ind2sub([n m], find(isfinite(A_j)));
        
    for k = 1:numel(idx_i)
        i = idx_i(k);
    
        if ~visited(i)
            visited(i) = 1;

            eta(i) = eta(j);
            v(i) = A_pi(i,j,idx_l(k)) - ((idx_l(k) - 1) * eta(j)) + v(j);
        
            [eta, v, visited] = visit_all(A_pi, xi, eta, v, i, visited);
        end
    end
end

function eta_bar = avg_weight(xi, A_pi)
    eta_bar = 0;
    len = 0;

    for i = 1:length(xi(:,1))-1
        eta_bar = eta_bar + A_pi(xi(i,1), xi(i+1,1), xi(i,2));
        len = len + (xi(i,2) - 1);
    end
    
    eta_bar = eta_bar / len;
end

function [eta, v, j] = value_det(A_pi, v_old)
    n = size(A_pi,1);
    %m = size(A_pi,3);

    eta = -inf * ones(n,1);
    v = -inf * ones(n,1);

    j = 0;

    s = 1;
    visited = zeros(1, n);

    while ~isempty(s)
        xi = cycle(A_pi, s);
        eta_bar = avg_weight(xi, A_pi);

        if j == 0, j = xi(1); end
        
        eta(xi(1,1)) = eta_bar;
        v(xi(1,1)) = v_old(xi(1,1));
        
        visited(xi(1,1)) = 1;

        [eta, v, visited] = visit_all(A_pi, xi(1:end-1,:), eta, v, xi(1,1), visited);

        s = find(visited == 0, 1);
    end
end

function [eta, v, j] = improvement(A, pi, A_pi, v)
    
    tol = 1e-10;

    n = size(A_pi,1);
    m = size(A_pi,3);
    cont_improvement = true;

    while cont_improvement

        cont_improvement1 = true;

        while cont_improvement1
            [eta, v, j] = value_det(A_pi, v);

            I_1 = 0;
            for i=1:n
                A_i = squeeze(A(i,:,:))';
                eta_j = repmat(eta', m, 1);

                eta_j(~isfinite(A_i)) = -inf;

                [max_eta, idx] = max(eta_j(:));
                [max_l, max_j] = ind2sub([m n], idx);
    
                if eta(i) < max_eta - tol
                    I_1 = I_1 + 1; 
                    A_pi(i,pi(i,1),pi(i,2)) = -inf;
    
                    pi(i,1) = max_j;
                    pi(i,2) = max_l;
                    A_pi(i,max_j,max_l) = A(i,max_j,max_l);
                end
            end
    
            if I_1 == 0
                P_bar = (eta(:) == eta(:).') .* isfinite(A);
                %P_bar = logical(P_bar);
                cont_improvement1 = false;
            end    
        end

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
            
            if v(i) < max_v - tol
                I_2 = I_2 + 1; 
                A_pi(i,pi(i,1),pi(i,2)) = -inf;
    
                pi(i,1) = max_j;
                pi(i,2) = max_l;
                A_pi(i,max_j,max_l) = A(i,max_j,max_l);
            end
        end
    
        if I_2 == 0, cont_improvement = false; end     
    end
end


