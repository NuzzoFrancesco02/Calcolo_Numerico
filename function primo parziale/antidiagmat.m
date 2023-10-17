function A=antidiagmat(x) 
%%crea una matrice di zeri con il vettore in input sull'antidiagonale
n=size(x,1);
A=zeros(n);
for j=1:n
    A(n+1-j,j)=x(j);
end
end