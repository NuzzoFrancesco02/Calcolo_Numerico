
function [L,U,x] = thomas(A,b)
    n = length(b);
    if (size(A,1)~=n || size(A,2)~=n)
        error('Dimensioni sistema errate');
    end
    for i = 2:n
        if diag(A,i)~=0 | diag(A,-i)~=0
            error('La matrice non è tridiagonale');
        end
    end
    alfa(1) = A(1,1);
    e = diag(A,-1);
    c = diag(A,1);
    for i = 2 : n
        d(i-1) = e(i-1)/alfa(i-1);
        alfa(i) = A(i,i)-d(i-1)*c(i-1);
    end
    L = eye(n) + diag(d,-1);
    U = diag(alfa) + diag(c,1);
    % LU*x=b --> U*x=y --> L*y=b
    y = zeros(n,1);
    y(1) = b(1);
    for i = 2:n
        y(i) = b(i)-d(i-1)*y(i-1);
    end
    x = zeros(n,1);
    x(n)=y(n)/alfa(n);
    for i = n-1:-1:1
        x(i) = (y(i)-c(i)*x(i+1))/alfa(i);
    end
end



