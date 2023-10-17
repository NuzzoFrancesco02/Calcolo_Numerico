%% NUMERO DI INTERVALLI PER AVERE UNA TOLLERNAZA
function N = integr_stima(a,b,toll,fun,str)
x = a:1e-5:b;
ddf = Jac(Jac(fun));
if strcmp(str,'pm')
    N = sqrt(((b-a)^3)*max(abs(ddf(x)))/(24*toll));
elseif strcmp(str,'t')
    N = sqrt(((b-a)^3)*max(abs(ddf(x)))/(12*toll));
elseif strcmp(str,'s')
    ddddf = Jac(Jac(ddf));
    N = nthroot(((b-a)^5)*max(abs(ddddf(x)))/(16*180*toll),4);
end