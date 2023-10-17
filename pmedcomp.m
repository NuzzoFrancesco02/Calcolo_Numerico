function I = pmedcomp(a,b,N,fun)
    H = (b-a)/N;
    x_nod = a+H/2:H:b-H/2;
    I = H*sum(fun(x_nod));
