function sdp(A)
if A==A'
    v=eig(A);
    if v>0
        disp('La matrice è SDP')
    else 
        disp('La matrice non è SDP')
    end
else disp('La matrice non è simmetrica')
end
end  