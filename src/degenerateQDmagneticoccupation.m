%Calculate the field due to the local QD magnetization for non-polarized QD
g0ret = 1./(w-epsilon+1i*g0tot/4);
g0adv = 1./(w-epsilon-1i*g0tot/4);
g0less = +1i*g0tot.*(fermi(1,:)+fermi(2,:))./((w-epsilon).^2+(g0tot/4)^2);

mx = -J/(2*pi).*imag(trapz(w, g0ret.*Sx(j).*g0less + g0less.*Sx(j).*g0adv));
my = -J/(2*pi).*imag(trapz(w, g0ret.*Sy(j).*g0less + g0less.*Sy(j).*g0adv));
mz = -J/(2*pi).*imag(trapz(w, g0ret.*Sz(j).*g0less + g0less.*Sz(j).*g0adv));
