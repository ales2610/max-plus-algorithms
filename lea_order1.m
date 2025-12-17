function [delta_min, delta_max, low_95, upp_95] = lea_order1(A, c, j, N_max, M, delta, f_sample, pr)

    n = size(A,1);

    delta_min = 0;
    delta_max = 0;

    avg_d = zeros(2,N_max);
    d_all = zeros(2,M);

    for m=1:M
        delta_m_min = -inf;
        delta_m_max = inf;

        A_k0 = identity(n);             
        A_kc = identity(n);  

        m

        for k=0:N_max-1
          
            A_k = f_sample(A);          

            A_k0 = otimes(A_k, A_k0);

            %k
            
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

    idx = find(isinf(avg_d(1, :)), 1, 'first');

    if pr
        x = idx:1:N_max;

        figure
        plot(x, avg_d(1, idx:end), x, avg_d(2, idx:end))

        %figure
        %histogram(d_all(1,1:end));
    end
end

function E = identity(n)
    E = -inf*ones(n);
    for i=1:n
        E(i,i) = 0;
    end
end

% function C = otimes(A, B)
%     n = size(A,1);
%     C = -inf*ones(n);
%     for i = 1:n
%         for j = 1:n
%             for k = 1:n
%                 C(i,j) = max(C(i,j), A(i,k) + B(k,j));
%             end
%         end
%     end
% end

function C = otimes(A, B)
    n = size(A,1);
    C = -inf(n);

    for k = 1:n
        C = max(C, A(:,k) + B(k,:));
    end
end