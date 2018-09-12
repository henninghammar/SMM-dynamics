function [G0ret, G0adv, G1xret, G1xadv, G1yret, G1yadv, G1zret, G1zadv] = stationaryretardedgreensfunction(pL, pR, gamma, mu, eps, w, J, S, wL, beta)

Sx = S(1);
Sy = S(2);
Sz = S(3);

g0L=gamma; %Gamma left
g0R=gamma; %Gamma righ
gSL=pL*g0L;
gSR=pR*g0R;
g(1)=(g0L/2*(1+pL)+g0R/2*(1+pR)); %Gamma up
g(2)=(g0L/2*(1-pL)+g0R/2*(1-pR)); %Gamma down

fermiL = 1./(1+exp(beta(1).*(w+mu(1))));
fermiR = 1./(1+exp(beta(2).*(w+mu(2))));

g0ret = 0.5*(1./(w-eps(1)+1i.*g(1)/2)+1./(w-eps(2)+1i.*g(2)/2));
g0adv = 0.5*(1./(w-eps(1)-1i.*g(1)/2)+1./(w-eps(2)-1i.*g(2)/2));
g1ret = 0.5*(1./(w-eps(1)+1i.*g(1)/2)-1./(w-eps(2)+1i.*g(2)/2));
g1adv = 0.5*(1./(w-eps(1)-1i.*g(1)/2)-1./(w-eps(2)-1i.*g(2)/2));

G0ret = g0ret - J*Sz.*(g0ret.*g1ret+g1ret.*g0ret);
G0adv = g0adv - J*Sz.*(g0adv.*g1adv+g1adv.*g0adv);
G1xret = - J.*(Sx.*g0ret.*g0ret-1i*Sy.*g1ret.*g0ret+1i*Sy.*g0ret.*g1ret+1i.*Sx.*g1ret.*g1ret);
G1xadv = - J.*(Sx.*g0adv.*g0adv-1i*Sy.*g1adv.*g0adv+1i*Sy.*g0adv.*g1adv+1i.*Sx.*g1adv.*g1adv);
G1yret = - J.*(Sy.*g0ret.*g0ret+1i*Sx.*g1ret.*g0ret-1i*Sx.*g0ret.*g1ret+1i.*Sy.*g1ret.*g1ret);
G1yadv = - J.*(Sy.*g0adv.*g0adv+1i*Sx.*g1adv.*g0adv-1i*Sx.*g0adv.*g1adv+1i.*Sy.*g1adv.*g1adv);
G1zret = g1ret - J*Sz.*g0ret.*g0ret;
G1zadv = g1adv - J*Sz.*g0adv.*g0adv;

end
