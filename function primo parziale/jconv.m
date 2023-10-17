function rho=jconv(A)

n=size(A,1);
D=diag(diag(A)); Dinv=diag(1./diag(A));
B=eye(n)-Dinv*A;
rho=spectre(B);
if spectre(B)<1 
    disp('Jacobi converge')
    c=1;
else
    disp('Jacobi non converge')
    c=0;
end
