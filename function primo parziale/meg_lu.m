function [x,L,U]=meg_lu(A,b)

n=size(A,1);
for k=1:n-1
    for i=k+1:n
        A(i,k)=A(i,k)/A(k,k);
        for j=k+1:n
            A(i,j)=A(i,j)-A(i,k)*A(k,j);
        end
        b(i)=b(i)-A(i,k)*b(k);
    end
    A
end

U=triu(A);
L=tril(A,-1)+eye(4,4);
x=bksub(U,b);