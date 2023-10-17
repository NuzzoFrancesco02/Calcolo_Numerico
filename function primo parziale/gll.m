function [nod,pes]=gll(n,a,b)
    [M,~]=zplege(n-2,a,b);
    nod=[-1;M;1];
    
    L=[];
    syms x real
    L1=1;
    L2=x;
    Lprecprec=L1;
    Lprec=L2;
    for i=3:n+2
        k=i-1;
        L= ((2*k+1)/(k+1))*x*Lprec-(k/(k+1))*Lprecprec;
        Lprecprec=Lprec;
        Lprec=L;
    end
    pes=[];
     x=nod;
     l=eval(L);
    pes=(2/(n*(n+1))).*(1/(l.^2));
   
end