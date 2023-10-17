function fn=es41(t,y)
 [n,m]=size(y);
    fn=zeros(n,m);
    fn(1)= -2*y(2) - 10*(y(1))^2 + exp(-t/2)*(2*cos(t)-(7/2)*sin(t))+40*exp(-t)*(sin(t))^2;
    fn(2)= y(2);
return
