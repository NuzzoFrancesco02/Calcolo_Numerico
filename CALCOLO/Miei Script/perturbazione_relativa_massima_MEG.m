clear
clc

A = [];
n = length(A);
epsm = 2^-52;

%%% MEG seza pivoting
maxK = max(max(A));
L = eye(n);
for k = 1:(n-1)				% fissa l'indice di colonna
    for i = (k+1):n			% fissa l'indice di riga       
        L(i,k) = A(i,k) / A(k,k);             
        for j = (k+1):n		% a partire dall'elemento fissato scorre la riga           
            A(i,j) = A(i,j) - L(i,k) * A(k,j);   
			maxK = max(max(max(abs(A))),maxK);
        end
    end
end
U = triu(A);
%%%

rhon = maxK/max(max(A));
pert = 8*(n^3)*epsm*rhon