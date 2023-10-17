function [xk, k, p, alpha] = conjgrad(A, b, x0, toll, nmax)
    x = x0;
    r = b - A*x;
    p = r;
    err = norm(r)/norm(b) + 1;
    [xk] = x0;
    k = 0;
    while err>toll && k<nmax
        k = k + 1;
        alpha = (p'*r)/(p'*A*p);
        x = x + alpha*p;
        r = r - alpha*A*p;
        beta = (p'*A*r)/(p'*A*p);
        p = r - beta * p;
        xk = [xk x];    
        err = norm(r)/norm(b);
    end
    
end