%% [D, k, err] = qrbasic (A, tol ,nmax)
% Autovalori non devono coincidere
% 
function [D, k, err] = qrbasic (A, tol ,nmax)
    if nargin == 1
        tol = 1e-6;
        nmax = 1000;
    end
    k = 0;
    err = 10;
    while err > tol && k < nmax
        k = k + 1;
        [Q, R] = qr(A);
        A = R * Q;
        D = diag(A);
        E = tril(A,-1);
        err = max(max(abs(E)));
    end
    if k >= nmax 
        fprintf('"qrbasic" non converge in %d iterazioni', k);
    else
        fprintf('"qrbasic" converge in %d iterazioni', k);
    end
end