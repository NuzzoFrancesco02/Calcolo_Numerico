function I = trapcomp(a,b,N,fun)
    H = (b-a)/N;
    xk = a+H:H:b-H;
    I = H*(fun(a)+fun(b))/2 + H*sum(fun(xk));