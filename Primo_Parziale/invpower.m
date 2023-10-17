function [lambda,x,iter] = invpower(A,x0,tol ,nmax)
    if nargin==1
        tol = 1e-6;
        x0 = ones(size(A,1),1);
        nmax = 1000;
    end
    if length(x0) ~= size(A,1) || length(x0) ~= size(A,2)
        error('Dimensioni errate')
    end
    
    iter = 0;
    y = x0/norm(x0);
    lambda = y' * A * y;
    [L, U, P] = lu(A);
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