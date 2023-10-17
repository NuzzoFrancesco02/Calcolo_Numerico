%% [lambda,x,iter,lambda_vect] = eigpower(A,tol ,nmax,x0,l_max)
% Autovalori diversi e diversi in modulo
function [lambda,x,iter,lambda_vect] = eigpower(A,tol ,nmax,x0)
    
    if length(x0) ~= size(A,1) || length(x0) ~= size(A,2)
        error('Dimensioni errate')
    end
    if nargin==1
        tol = 1e-6;
        x0 = ones(n,1);
        nmax = 1000;
    end
    
    matr = [];
    iter = 0;
    y = x0/norm(x0);
    lambda = y' * A * y;
    lambda_vect = [lambda];
    
    err = tol * abs(lambda) + 1;
    while err > tol * abs(lambda) && abs(lambda) ~= 0 && iter < nmax
        iter = iter + 1;
        x = A * y;
        y = x/norm(x);
        matr = [matr, y];
        lambdanew = y' * A * y;
        err = norm(lambdanew-lambda);
        lambda = lambdanew;
        lambda_vect = [lambda_vect lambda];
    end
    if (err <= tol * abs(lambda))
        fprintf('eigpower converge in %d iterazioni \n', iter);
    else
        fprintf('eigpower non converge in %d iterazioni \n', iter);
    end
end