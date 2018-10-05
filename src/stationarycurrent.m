function [Ic, Isx, Isy, Isz, IeL, IqL, IeR, IqR, InL, InR, DOS, MDOSz, Iq0L, Iq1L, Iq0R, Iq1R, Ic0, Ic1, Is0z, Is1z, Ie0L, Ie1L] = stationarycurrent(pL, pR, gamma, mu, eps, w, dw, J, S, wL, beta)
    [G0less, G0great, G1xless, G1xgreat, G1yless, G1ygreat, G1zless, G1zgreat] = stationarygreensfunction(pL, pR, gamma, mu, eps, w, J, S, wL, beta);

    g0L=gamma; %Gamma left
    g0R=gamma; %Gamma righ
    gSL=pL*g0L;
    gSR=pR*g0R;
    g(1)=(g0L/2*(1+pL)+g0R/2*(1+pR)); %Gamma up
    g(2)=(g0L/2*(1-pL)+g0R/2*(1-pR)); %Gamma down

    fermiL = 1./(1+exp(beta(1).*(w+mu(1))));
    fermiR = 1./(1+exp(beta(2).*(w+mu(2))));

    self0less = 1i*g0L*fermiL;
    self0great = -1i*g0L*(1-fermiL);
    self1less = 1i*gSL*fermiL;
    self1great = -1i*gSL*(1-fermiL);

    self0lessR = 1i*g0R*fermiR;
    self0greatR = -1i*g0R*(1-fermiR);
    self1lessR = 1i*gSR*fermiR;
    self1greatR = -1i*gSR*(1-fermiR);

    %Charge current
    Ic0bare = -(1/pi)*(self0less.*G0great - self0great.*G0less);
    Ic1bare = -(1/pi)*(self1less.*G1zgreat - self1great.*G1zless);
    Ic0 = trapz(w,Ic0bare);
    Ic1 = trapz(w,Ic1bare);
    Ic = Ic0 + Ic1;

    %Spin current
    Is0xbare = -(1/pi)*(self0less.*G1xgreat - self0great.*G1xless);
    Is0ybare = -(1/pi)*(self0less.*G1ygreat - self0great.*G1yless);
    Is0zbare = -(1/pi)*(self0less.*G1zgreat - self0great.*G1zless);
    Is1xbare = 1i*(1/pi)*(self1less.*G1ygreat - self1great.*G1yless);
    Is1ybare = -1i*(1/pi)*(self1less.*G1xgreat - self1great.*G1xless);
    Is1zbare = -(1/pi)*(self1less.*G0great - self1great.*G0less);
    Is0x = trapz(w,Is0xbare);
    Is0y = trapz(w,Is0ybare);
    Is0z = trapz(w,Is0zbare);
    Is1x = trapz(w,Is1xbare);
    Is1y = trapz(w,Is1ybare);
    Is1z = trapz(w,Is1zbare);
    Isx = Is0x + Is1x;
    Isy = Is0y + Is1y;
    Isz = Is0z + Is1z;

    %Particle current
    In0bare = (1/pi)*(self0less.*G0great - self0great.*G0less);
    In1bare = (1/pi)*(self1less.*G1zgreat - self1great.*G1zless);
    In0L = trapz(w,In0bare);
    In1L = trapz(w,In1bare);
    InL = In0L + In1L;

    In0bareR = (1/pi)*(self0lessR.*G0great - self0greatR.*G0less);
    In1bareR = (1/pi)*(self1lessR.*G1zgreat - self1greatR.*G1zless);
    In0R = trapz(w,In0bareR);
    In1R = trapz(w,In1bareR);
    InR = In0R + In1R;

    %Energy current
    Ie0bare = (1/pi).*(w- mu(1)).*(self0less.*G0great - self0great.*G0less);
    Ie1bare = (1/pi).*(w- mu(1)).*(self1less.*G1zgreat - self1great.*G1zless);
    Ie0L = trapz(w,Ie0bare);
    Ie1L = trapz(w,Ie1bare);
    IeL = Ie0L + Ie1L;

    Ie0bareR = (1/pi).*(w- mu(2)).*(self0lessR.*G0great - self0greatR.*G0less);
    Ie1bareR = (1/pi).*(w- mu(2)).*(self1lessR.*G1zgreat - self1greatR.*G1zless);
    Ie0R = trapz(w,Ie0bareR);
    Ie1R = trapz(w,Ie1bareR);
    IeR = Ie0R + Ie1R;

    %Heat current
    Iq0L = Ie0L - mu(1)*In0L;
    Iq0R = Ie0R - mu(2)*In0R;
    Iq1L = Ie1L - mu(1)*In1L;
    Iq1R = Ie1R - mu(2)*In1R;
    IqL = Iq0L + Iq1L;
    IqR = Iq0R + Iq1R;

    %In SI units
    Iconv=1.602176565*10^(-19)/(6.58211928*10^(-16))*10^-3;%in A, elementary charge divided by hbar in eVs times 10^-3 as we use meV
    Ic = Ic.*Iconv;
    Isx = Isx.*Iconv;
    Isy = Isy.*Iconv;
    Isz = Isz.*Iconv;

    %energyconv = 1/(6.58211928*10^(-16))*10^-3;%1 divided by hbar in eVs times 10^-3 as we use meV
    %Ie = Ie*energyconv;
    %Iq = Iq*energyconv;

    DOS = 1i/(2*pi)*(G0great-G0less);
    MDOSz = 1i/(4*pi)*(G1zgreat-G1zless);
end
