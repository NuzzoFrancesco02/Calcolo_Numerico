function MU=mu_upwind_ordine2(mu,eta,h)
  Pe=pechlet(mu,eta,h);
  MU=mu*(1+(Pe-1+((2*Pe)/(exp(2*Pe))-1)));
end
