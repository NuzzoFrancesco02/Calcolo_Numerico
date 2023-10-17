function [L,U] = lugauss(A)

% function [L,U] = lugauss(A)
% 
% A: matrice quadrata da fattorizzare
% la matrice L ha tutti 1 sulla diagonale

[n,m] = size(A);

if (n ~= m)
    error('A non è una matrice quadrata'); 
end

L = eye(n);

for k = 1:(n-1) % fissa l'indice di colonna
    
    if (A(k,k) == 0)
        error('Un elemento pivot si è annullato'); 
    end
    
    for i = (k+1):n % fissa l'indice di riga               
        
        L(i,k) = A(i,k) / A(k,k);       
        
        for j = (k+1):n % a partire dall'elemento fissato scorre la riga
            
            A(i,j) = A(i,j) - L(i,k) * A(k,j);            
        end
    end
end

U = triu(A);
end