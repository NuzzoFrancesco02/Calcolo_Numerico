function [lambda,x,iter] = invpowershift(A, s, tol ,nmax,x0)
    
    if nargin==2
        tol = 1e-6;
        x0 = ones(size(A,1),1);
        nmax = 1000;
    end
    n = length(x0);
    if n ~= size(A,1) || n ~= size(A,2)
        error('Dimensioni errate')
    end
    
    iter = 0;
    y = x0/norm(x0);
    lambda = y' * A * y;
    M = A - s*eye(n);
    [L, U, P] = lu(M);
    err = tol * abs(lambda) + 1;
    while err > tol * abs(lambda) && abs(lambda) ~= 0 && iter < nmax
        iter = iter + 1;
        c = L \ P*y;
        x = U \ c;
        y = x/norm(x);
        lambdanew = y' * A * y;
        err = norm(lambdanew-lambda);
        lambda = lambdanew;
    end
    if (err <= tol * abs(lambda))
        fprintf('invpower converge in %d iterazioni \n', iter);
    else
        fprintf('invpower non converge in %d iterazioni \n', iter);
    end
end