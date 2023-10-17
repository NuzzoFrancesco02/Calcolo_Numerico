function MU=mu_upwind_ordine1(mu,eta,h)
    MU=mu*(1+pechlet(mu,eta,h));
end