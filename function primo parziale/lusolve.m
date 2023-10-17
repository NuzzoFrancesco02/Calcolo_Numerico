function X=lusolve(A,B)
%risolve sistemi lineari con la stessa matrice dei coefficienti

[L,U]=lufat(A);
n=size(A,2);
s=size(B,1);
Y=ones(s,n); X=ones(s,n);

for i=1:s
    Y(i,:)=(fwsub(L,B(i,:)'))';
end

for i=1:s
    X(i,:)=(bksub(U,Y(i,:)'))';
end
end