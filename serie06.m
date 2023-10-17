%% serie 06
V = pi*((sinh(1)+sinh(3/5))/4 + 2/5);
f = @(x) pi*cosh(x-0.5).^2;
a = 0;
b = 0.8;
x = 0:0.001:0.8;
for N = 1:20
    I(1,N) = pmedcomp(0,0.8,N,f);
    I(2,N) = trapcomp(0,0.8,N,f);
    I(3,N) = simpcomp(0,0.8,N,f);
    E(1,N) = abs(V-I(1,N));
    E(2,N) = abs(V-I(2,N));
    E(3,N) = abs(V-I(3,N));
    E_gauss(N) = abs(V-gausscomp(0,0.8,N,f));
    H(N) = (b-a)/N;
end
N = 1:20;

stimap_2(E(1,:),H);
figure()
plot(1:20,I,'o','MarkerSize',5,'LineWidth',2)
legend('Pto medio','Trapez','Simpson')
hold on
plot(1:20,V*ones(20,1),'LineWidth',2)
figure()
loglog(H,E,'LineWidth',2);
hold on;
loglog(H,H,'-k',H,H.^2,'-.k',H,H.^3,':k',H,H.^4,'k')
legend('Pto medio','Trapez','Simpson','H','H^2','H^3','H^4')

N_pm = ceil(integr_stima(a,b,1e-5,f,'pm'))
abs(V-pmedcomp(a,b,N_pm,f))
N_t = ceil(integr_stima(a,b,1e-5,f,'t'))
abs(V-trapcomp(a,b,N_t,f))
N_s = ceil(integr_stima(a,b,1e-5,f,'s'))
abs(V-simpcomp(a,b,N_s,f))

figure()
loglog(H,E(2,:),'r',H,E(3,:),'b',H,E_gauss,'y','LineWidth',2)
hold on
loglog(H,H,'-k',H,H.^2,'-.k',H,H.^4,':k')
legend('Trapezio','Simpson','Gauss composito','H','H^2','H^4')
p = stimap_2(E(3,:),H);
p(end)
%% 3
p1 = @(x) x.^4 -2*x + 1;
p2 = @(x) 3*x.^9 -5*x.^4 + 1;
p3 = @(x) 10./(x+2);
p4 = @(x) sqrt(x);

nv = 0:7;
Ep1_gauss = [];
for n = nv
    [xi,wi] = zplege(n,0,1);
    I = sum(wi.*p1(xi));
    Ep1_gauss = [Ep1_gauss, abs(I-0.2)];
end
Ep2_gauss = [];
for n = nv
    [xi,wi] = zplege(n,0,1);
    I = sum(wi.*p1(xi));
    Ep2_gauss = [Ep2_gauss, abs(I-0.3)];
end
Ep3_gauss = [];
for n = nv
    [xi,wi] = zplege(n,0,1);
    I = sum(wi.*p1(xi));
    Ep3_gauss = [Ep3_gauss, abs(I-10*log(3/2))];
end
Ep4_gauss = [];
for n = nv
    [xi,wi] = zplege(n,0,1);
    I = sum(wi.*p1(xi));
    Ep4_gauss = [Ep4_gauss, abs(I-2/3)];
end



semilogy(nv,max(Ep1_gauss, 1e-15),'LineWidth',2)
hold on;
semilogy(nv,max(Ep2_gauss, 1e-15),'LineWidth',2)
semilogy(nv,max(Ep3_gauss, 1e-15),'LineWidth',2)
semilogy(nv,max(Ep4_gauss, 1e-15),'LineWidth',2)



