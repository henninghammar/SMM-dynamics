Sx = S(1);
Sy = S(2);
Sz = S(3);

g0L=gamma; %Gamma left
g0R=gamma; %Gamma righ
gSL=pL*g0L;
gSR=pR*g0R;
g(1)=(g0L/2*(1+pL)+g0R/2*(1+pR)); %Gamma up
g(2)=(g0L/2*(1-pL)+g0R/2*(1-pR)); %Gamma down

fermiL = 1./(1+exp(beta(1).*(w+eV(1))));
fermiR = 1./(1+exp(beta(2).*(w+eV(2))));

g0ret = 0.5*(1./(w-eps(1)+1i.*g(1)/2)+1./(w-eps(2)+1i.*g(2)/2));
g0adv = 0.5*(1./(w-eps(1)-1i.*g(1)/2)+1./(w-eps(2)-1i.*g(2)/2));
g1ret = 0.5*(1./(w-eps(1)+1i.*g(1)/2)-1./(w-eps(2)+1i.*g(2)/2));
g1adv = 0.5*(1./(w-eps(1)-1i.*g(1)/2)-1./(w-eps(2)-1i.*g(2)/2));

g0less = 0.5*1i*((g0L*fermiL+g0R*fermiR).*(1./((w-eps(1)).^2+(g(1)/2)^2)+1./((w-eps(2)).^2+(g(2)/2)^2))...
  +(gSL*fermiL+gSR*fermiR).*(1./((w-eps(1)).^2+(g(1)/2)^2)-1./((w-eps(2)).^2+(g(2)/2)^2)));
g1less = 0.5*1i*((gSL*fermiL+gSR*fermiR).*(1./((w-eps(1)).^2+(g(1)/2)^2)+1./((w-eps(2)).^2+(g(2)/2)^2))...
  +(g0L*fermiL+g0R*fermiR).*(1./((w-eps(1)).^2+(g(1)/2)^2)-1./((w-eps(2)).^2+(g(2)/2)^2)));
g0great = -0.5*1i*((g0L*(1-fermiL)+g0R*(1-fermiR)).*(1./((w-eps(1)).^2+(g(1)/2)^2)+1./((w-eps(2)).^2+(g(2)/2)^2))...
  +(gSL*(1-fermiL)+gSR*(1-fermiR)).*(1./((w-eps(1)).^2+(g(1)/2)^2)-1./((w-eps(2)).^2+(g(2)/2)^2)));
g1great = -0.5*1i*((gSL*(1-fermiL)+gSR*(1-fermiR)).*(1./((w-eps(1)).^2+(g(1)/2)^2)+1./((w-eps(2)).^2+(g(2)/2)^2))...
  +(g0L*(1-fermiL)+g0R*(1-fermiR)).*(1./((w-eps(1)).^2+(g(1)/2)^2)-1./((w-eps(2)).^2+(g(2)/2)^2)));

G0less = g0less - J*Sz.*(g0ret.*g1less+g0less.*g1adv+g1ret.*g0less+g1less.*g0adv);
G0great = g0great - J*Sz.*(g0ret.*g1great+g0great.*g1adv+g1ret.*g0great+g1great.*g0adv);
G1xless = - J.*(Sx.*(g0ret.*g0less+g0less.*g0adv)-1i*Sy.*(g1ret.*g0less+g1less.*g0adv)...
  +1i*Sy.*(g0ret.*g1less+g0less.*g1adv)+1i.*Sx.*(g1ret.*g1less+g1less.*g1adv));
G1xgreat = - J.*(Sx.*(g0ret.*g0great+g0great.*g0adv)-1i*Sy.*(g1ret.*g0great+g1great.*g0adv)...
  +1i*Sy.*(g0ret.*g1great+g0great.*g1adv)+1i.*Sx.*(g1ret.*g1great+g1great.*g1adv));
G1yless = - J.*(Sy.*(g0ret.*g0less+g0less.*g0adv)+1i*Sx.*(g1ret.*g0less+g1less.*g0adv)...
  -1i*Sx.*(g0ret.*g1less+g0less.*g1adv)+1i.*Sy.*(g1ret.*g1less+g1less.*g1adv));
G1ygreat = - J.*(Sy.*(g0ret.*g0great+g0great.*g0adv)+1i*Sx.*(g1ret.*g0great+g1great.*g0adv)...
  -1i*Sx.*(g0ret.*g1great+g0great.*g1adv)+1i.*Sy.*(g1ret.*g1great+g1great.*g1adv));
G1zless = g1less - J*Sz.*(g0ret.*g0less+g0less.*g0adv);
G1zgreat = g1great - J*Sz.*(g0ret.*g0great+g0great.*g0adv);
