% Controllo la condizione necessaria e sufficiente
function boolean = cond_nec_suff(M)
    boolean = false;
    n = length(diag(M));
    det_m = 1;
    i = 1;
    while ( i < n && det_m ~= 0)
        det_m = det(M(1:i,1:i));
        i = i +1;
    end
    if (det_m ~= 0)
        boolean = true;
    else
        error('Non è soddisfatta la condizione necessaria e sufficiente!');
    end
end
