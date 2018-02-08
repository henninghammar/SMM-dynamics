function [mx, my, mz] = degenerateQDmagneticoccupation(t, t0, t1, epsilon, g0tot, g0, mu, w, fermi, S, J)

[G0less, G0great, G1xless, G1xgreat, G1yless, G1ygreat, G1zless, G1zgreat] = SMMnonpolarizedgreensfunction(t, t, t0, t1, epsilon, g0tot, g0, mu, w, fermi, S, J);

mx = imag(G1xless);
my = imag(G1yless);
mz = imag(G1zless);

end
