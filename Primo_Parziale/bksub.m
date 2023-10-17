%[y] = bksub(U,x)
function [x] = fwsub(U, b)

    n = length(b);
    if ((size(U,1) ~= n)|(size(U,2) ~= n))
        error ('Dimensioni non compatibili');
    end
    if ~isequal(U,triu(U))
        error('La matrice non è triangolare superiore');
    end
    if (prod(diag(U))==0)
        error('Almeno un elemento della diagonale è nullo');
    end
    x = zeros(n,1);
    x(n) = b(n)/U(n,n);
    for i = n-1 :-1: 1
        x(i) = (b(i) - U(i,i+1:n)*x(i+1:n))/U(i,i);
    end
end