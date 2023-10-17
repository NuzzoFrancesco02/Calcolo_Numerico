function [x,k] = jacobiconinv(A,b,x0,toll,nmax)
n = size(A,1);
D = diag(A);
e = diag(eye(n));
D_inversa = diag(e./D);
x = x0;
r = b - A*x;
k = 1;
while norm(r)./norm(b) > toll && k < nmax
    z = D_inversa * r; % Non bisogna calcolare inv perche` per una diagonale basta dividere elementi su diagonale per uno
    x = x +z;
    k = k+1;
    r = b - A*x;
end
end


    