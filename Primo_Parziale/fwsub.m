%[x] = fwsub(L,b)
function [x] = fwsub(L, b)
    
    n = length(b);
    if ((size(L,1) ~= n)|(size(L,2)~= n))
        error ('Dimensioni non compatibili');
    end
    if ~isequal(L,tril(L))
        error('La matrice non è triangolare inferiore');
    end
    if (prod(diag(L))==0)
        error('Almeno un elemento della diagonale è nullo');
    end
    x = zeros(n,1);
    
    x(1) = b(1)/L(1,1);
    for i = 2 : n
        x(i) = (b(i) - L(i,1:i-1)*x(1:i-1))/L(i,i);
    end
end
