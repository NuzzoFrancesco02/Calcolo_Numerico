function c=richarconv(A,P,alpha)

lambda=(eig(P\A));
c=[];
for i=1:length(lambda)
    if 2*real(lambda(i))/(alpha*abs(lambda(i)^2))>1
        c=[c 1];
    else
        c=[c 0];
    end
end

if c==ones(1,length(lambda))
    disp('Richardson stazionario converge')
else
    disp('Richardson stazionario non converge')
end
end