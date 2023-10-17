function [L,U] = lugauss(A)
    n = size(A,1);
    if ~isequal(size(A,1),size(A,2))
        error('Matrice non quadrata');
    end
    L = eye(n);
    for k = 1 : n - 1 
        if A(k,k) == 0
            error('Un elemento pivotale si è annullato');
        end
        for i = k + 1 : n
            L(i,k) = A(i,k)/A(k,k);
             for j = k +1 : n
                 A(i,j) = A(i,j) - L(i,k)*A(k,j);
             end
            %A(i,k:n) = A(i,k:n)-A(k,k:n)*L(i,k);
        end
   
    end
    U = triu(A);
end