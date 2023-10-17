function lu_necsuf(A)

i=1;
da=1;
na=size(A,1);

while ((i<na) && (da~=0))
    da=det(A(1:i,1:i));
    i=i+1;
end

if da==0
    disp('LU non applicabile')
else
    disp('LU applicabile')
end
end