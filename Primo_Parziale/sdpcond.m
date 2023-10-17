function K = sdpcond(A, tol ,nmax)
n = size(A,1);
if ((A ~= A') | (min(eig(A))<=0))
    error('La matrice non è simmetrica e def. positiva')
end
l_max = eigpower(A,tol,nmax,ones(n,1));
l_min = invpower(A,tol,nmax,ones(n,1));
K = l_max/l_min;
