function [L,U]=lufat(A)
%fattorizza A in due matrici L e U
n=size(A,1);
for k=1:n-1
    for i=k+1:n
        A(i,k)=A(i,k)/A(k,k);
        for j=k+1:n
            A(i,j)=A(i,j)-A(i,k)*A(k,j);
        end
    end
end

L=tril(A,-1)+eye(n);
U=triu(A);
end