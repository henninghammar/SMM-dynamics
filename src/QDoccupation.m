function n = QDmagneticoccupation(t, t0, t1, eps, g, g0, gS, mu, w, fermi, S, J)

[G0less, G0great, G1xless, G1xgreat, G1yless, G1ygreat, G1zless, G1zgreat] = greensfunction(t, t, t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

n = 2*imag(G0less);

end
