function [delta_min, delta_max, low_95, upp_95] = lea_order2(A0, A1, c, j, N_max, M, delta, f_sample, pr)

    n = size(A0,1);

    delta_min = 0;
    delta_max = 0;

    avg_d = zeros(2,N_max);

    d_all = zeros(2,M);

    for m=1:M
        delta_m_min = -inf;
        delta_m_max = inf;

        A_k0 = identity(n);             
        A_kc = identity(n);  

        for k=0:N_max-1
          
            A1_k = f_sample(A1);          
            A0_k = f_sample(A0);     
            A_k = one_A(A0_k, A1_k);
        
            A_k0 = otimes(A_k, A_k0);
            
            if k >= c
                A_kc = otimes(A_k, A_kc);
            end

            Delta = A_k0 - A_kc;

            delta_m_min = min(Delta(:,j));
            if isnan(delta_m_min), delta_m_min = -inf; end
     
            delta_m_max = max(Delta(:,j));
            if isnan(delta_m_max), delta_m_max = inf; end

            avg_d(1, k+1) = avg_d(1, k+1) + delta_m_min;
            avg_d(2, k+1) = avg_d(2, k+1) + delta_m_max;

            if ~pr && delta_m_max - delta_m_min <= delta, break; end
        end

        delta_min = delta_min + delta_m_min;
        delta_max = delta_max + delta_m_max;

        d_all(1,m) = delta_m_min;
        d_all(2,m) = delta_m_max;
    end

    delta_min = delta_min / (c * M);
    delta_max = delta_max / (c * M);

    avg_d = avg_d ./ (c * M);

    d_all = d_all ./ c;

    low_95 = prctile(d_all(1,:), 2.5);
    upp_95 = prctile(d_all(1,:), 97.5);

    if pr
        x = 1:1:N_max;

        figure
        plot(x, avg_d(1, 1:end), x, avg_d(2, 1:end))

        histogram(d_all(1,1:end));
    end
end

function A = one_A(A0, A1)

    n = size(A0,1);
    A0_star = identity(n);
    A0_max = A0_star;

    for l=1:n-1
        A0_star = otimes(A0_star, A0);
        A0_max = max(A0_max, A0_star);
    end

    A = otimes(A0_max, A1);
end

function E = identity(n)
    E = -inf*ones(n);
    for i=1:n
        E(i,i) = 0;
    end
end

function C = otimes(A, B)
    n = size(A,1);
    C = -inf(n);

    for k = 1:n
        C = max(C, A(:,k) + B(k,:));
    end
end
