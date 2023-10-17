%{
function succ = treuno(n)
succ(1,1) = n;
i = 1;
while((succ(1,i) ~= 1))
    if(mod(succ(1,i),2) == 0)
        succ(1,i+1) = succ(1,i)/2;
    else
        succ(1,i+1) = 3.*succ(1,i) + 1;
    end
    i=i+1;
end
return
%}
function succ = treuno(n)
succ = [ n ];

while((succ(end) ~= 1))
    if(mod(succ(end),2) == 0)
        succ(end + 1) = succ(end)/2;
    else
        succ(end + 1) = 3.*succ(end) + 1;
    end
    
end
return