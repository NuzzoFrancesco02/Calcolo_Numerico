%% Matrice multi-diagonale 
%
%   M = diagonals(v,n)
%
%   # INPUT
%   v : vettore che contiene i coefficienti delle diagonali
%   n : dimensione matrice quadrata
function M = diagonals(v,n)
    u = length(v);
    M = zeros(n);
    if mod(u,2)~=0
        iter = ceil(u/2);
        ind1 = iter;
        ind2 = iter;
        for i = 0 : iter-1
            if i == 0
                M = M + diag(v(ind1)*ones(n-i,1),i);
            else
                M = M + diag(v(ind1)*ones(n-i,1),i);
                M = M + diag(v(ind2)*ones(n-i,1),-i);
            end
            ind1 = ind1 + 1;
            ind2 = ind2 -1;
        end
    else
        iter = ceil(u/2);
        v = [v(1:iter) 0 v(iter+1:end)];
        iter = iter + 1;
        ind1 = iter;
        ind2 = iter;
        for i = 0 : iter-1
            M = M + diag(v(ind1)*ones(n-i,1),i);
            M = M + diag(v(ind2)*ones(n-i,1),-i);
            ind1 = ind1 + 1;
            ind2 = ind2 -1;
        end
    end
end
