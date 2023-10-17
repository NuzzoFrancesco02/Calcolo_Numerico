function x=meg_smart(A,b)

n=size(A,1);
for k=1:n-1
    for i=k+1:n
        l_ik=A(i,k)/A(k,k);
        for j=k:n
            A(i,j)=A(i,j)-l_ik*A(k,j);
        end
        b(i)=b(i)-l_ik*b(k);
    end
    A
end

U=triu(A);
x=bksub(A,b);