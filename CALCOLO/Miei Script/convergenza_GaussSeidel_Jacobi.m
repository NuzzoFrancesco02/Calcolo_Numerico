syms x;
A = [];

D = diag(diag(A));
L = tril(A);
Bj = eye(size(A)) - D^-1*A;
Bgs = eye(size(A)) - L^-1*A;
rhoBj = max(abs(eig(Bj)))		% <1 converge sempre
rhoBgs = max(abs(eig(Bgs)))		% <1 converge sempre

%R_j = -log(rhoBj)
%R_gs = -log(rhoBgs)

%abs(eig(Bj))
%abs(eig(Bgs))

%gauss_seidel(A,b,x0,tol,nmax)
%jacobi(A,b,x0,tol,nmax)