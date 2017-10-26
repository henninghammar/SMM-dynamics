function [mx, my, mz] = constantQDmagneticoccupation(eps, gamma, gamma0, gammaS, eV, w, fermi, S, J)

%Calculate the field due to the local QD magnetization for polarized QD
g0ret = 0.5*(1./(w-eps(1)+1i*gamma(1)/2) + 1./(w-eps(1)+1i*gamma(1)/2));
g0adv = 0.5*(1./(w-eps(1)-1i*gamma(1)/2) + 1./(w-eps(1)-1i*gamma(1)/2));
g1ret = 0.5*(1./(w-eps(1)+1i*gamma(1)/2) - 1./(w-eps(1)+1i*gamma(1)/2));
g1adv = 0.5*(1./(w-eps(1)-1i*gamma(1)/2) - 1./(w-eps(1)-1i*gamma(1)/2));
g0less = +1i*0.5*((gamma0(1)*fermi(1,:)+gamma0(2)*fermi(2,:)).*(1./((w-eps(1)).^2+(gamma(1)/2)^2)+1./((w-eps(2)).^2+(gamma(2)/2)^2))...
+(gammaS(1)*fermi(1,:)+gammaS(2)*fermi(2,:)).*(1./((w-eps(1)).^2+(gamma(1)/2)^2)-1./((w-eps(2)).^2+(gamma(2)/2)^2)));
g1less = +1i*0.5*((gammaS(1)*fermi(1,:)+gammaS(2)*fermi(2,:)).*(1./((w-eps(1)).^2+(gamma(1)/2)^2)+1./((w-eps(2)).^2+(gamma(2)/2)^2))...
+(gamma0(1)*fermi(1,:)+gamma0(2)*fermi(2,:)).*(1./((w-eps(1)).^2+(gamma(1)/2)^2)-1./((w-eps(2)).^2+(gamma(2)/2)^2)));

G1xless = -J*(S(1)*(g0ret.*g0less + g0less.*g0adv) - 1i*S(2)*(g1ret.*g0less+g1less.*g0adv) + 1i*S(2)*(g0ret.*g1less+g0less.*g1adv) + 1i*S(1).*(g1ret.*g1less+g1less.*g1adv));
G1yless = -J*(S(2)*(g0ret.*g0less + g0less.*g0adv) + 1i*S(1)*(g1ret.*g0less+g1less.*g0adv) - 1i*S(1)*(g0ret.*g1less+g0less.*g1adv) + 1i*S(2).*(g1ret.*g1less+g1less.*g1adv));
G1zless = g1less - J*S(3)*(g0ret.*g0less + g0less.*g0adv+g1ret.*g1less+g1less.*g1adv);

mx = 1/(2*pi).*imag(trapz(w, G1xless));
my = 1/(2*pi).*imag(trapz(w, G1yless));
mz = 1/(2*pi).*imag(trapz(w, G1zless));

end
