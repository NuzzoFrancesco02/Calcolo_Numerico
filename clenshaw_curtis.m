function I = clenshaw_curtis(fun,n)
    if mod(n,2) == 0
        for k = 0 : n
            f = @(t) fun(cos(t)).*cos(k.*t);
            a(k+1) = 2.*(simpcomp(0,pi,1e5,f))./pi;
        end
        I = a(1);
        for k = 1 : n/2
            I = I + a(1+2*k)/(1-(2*k)^2);
        end
    else
        error('n non pari')
    end