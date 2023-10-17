%% [err_stim, err_vero] = perturbated_matrix(A,x,b,db)
function [err_stim, err_vero] = perturbated_matrix(A,x,b,db)
    % dA = 0 --> r = -db 
    err_stim = cond(A)*norm(-db)/norm(b);
    
    x_dx = A \ (b+db);
    err_vero = norm(x-x_dx)/norm(x);