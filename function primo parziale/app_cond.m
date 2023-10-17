function [condvett]=app_cond(A,tol,nmax,x0)
    [n,m] = size(A);
if n ~= m
    error('Solo per matrici quadrate');
end

if nargin == 1
   tol = 1.e-06;   x0 = ones(n,1);   nmax = 100;
end

% iterazione zero fuori dal ciclo while
iter = 0;
K=0;
c=cond(A);
condvett=[];
ymax= x0/norm(x0); % y0
ymin=ymax;
err = tol + 1; % dobbiamo entrare nel ciclo
while(err>tol && iter<nmax)
    iter=iter+1;
    xmax=A*ymax;
    ymax=xmax/norm(xmax);
    xmin=A\ymin;
    ymin=xmin/norm(xmin);
    K=(ymax'*A*ymax)/(ymin'*A*ymin);
    condvett=[condvett;K];
    err=abs(K-c);
end
if (err <= tol)
     fprintf(' converge in %d iterazioni \n', iter);
else
     fprintf(' non converge in %d iterazioni. \n', iter)
end

return
end








    



