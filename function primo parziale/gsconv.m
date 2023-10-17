function rho=gsconv(A)

n=size(A,1);
T=tril(A);
B=eye(n)-T\A;
rho=spectre(B);
if spectre(B)<1 
    disp('Gauss-Seidel converge')
else
    disp('Gauss-Seidel non converge')
end