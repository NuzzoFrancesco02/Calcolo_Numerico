function x=antidiag(A) 
%%estrae l'antidiagonale dalla matrice in input
n=size(A);
n=n(1);
x=zeros(n,1);
    for i=n:-1:1
        x(n+1-i)=A(i,n+1-i);
    end
end

