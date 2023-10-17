% function stim = CGL_stim(f,a,b,n)
%     A = lebesgue_CGL(n,a);
%     x = linspace(a,b,1e6);
%     stim = max(abs(f(x)-lagr_polin(n,f,a,b,x,'CGL')));