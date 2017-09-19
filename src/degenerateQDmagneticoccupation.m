function [mx, my, mz] = degenerateQDmagneticoccupation(w, epsilon, g0tot, J, fermi, S)

%Calculate the field due to the local QD magnetization for non-polarized QD
g0ret = 1./(w-epsilon+1i*g0tot/4);
g0adv = 1./(w-epsilon-1i*g0tot/4);
g0less = +1i*g0tot.*(fermi(1,:)+fermi(2,:))./((w-epsilon).^2+(g0tot/4)^2);

mx = -J/(2*pi).*imag(trapz(w, g0ret.*S(1).*g0less + g0less.*S(1).*g0adv));
my = -J/(2*pi).*imag(trapz(w, g0ret.*S(2).*g0less + g0less.*S(2).*g0adv));
mz = -J/(2*pi).*imag(trapz(w, g0ret.*S(3).*g0less + g0less.*S(3).*g0adv));

end
