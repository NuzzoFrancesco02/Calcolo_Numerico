clear all
clc

% esercizio 2

A = [ 50 1 3; 1 6 0; 3 0 1 ]
n = size(A,1);

[L,U] = lugauss(A);

I = eye(n);

Ainv = [];

for i = 1:n 
    yn = fwsub( L, I(:,i));
    xn = bksub( U, yn );
    Ainv = [Ainv, xn];
end
    
Ainv
eye(3) - Ainv * A