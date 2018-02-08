function [Ic, Isx, Isy, Isz, Ie, Iq, DOS, MDOSz] = stationarycurrent(pL, pR, gamma, eV, eps, epsilon, w, dw, J, S, wL, beta)
    stationarygreensfunction

    self0less = 1i*g0L*fermiL;
    self0great = -1i*g0L*(1-fermiL);
    self1less = 1i*gSL*fermiL;
    self1great = -1i*gSL*(1-fermiL);

    %Charge current
    Ic0bare = -(1/(2*pi))*(self0less.*G0great - self0great.*G0less);
    Ic1bare = -(1/(2*pi))*(self1less.*G1zgreat - self1great.*G1zless);
    Ic0 = trapz(w,Ic0bare);
    Ic1 = trapz(w,Ic1bare);
    Ic = Ic0 + Ic1;

    %Spin current
    Is0xbare = -(1/(2*pi))*(self0less.*G1xgreat - self0great.*G1xless);
    Is0ybare = -(1/(2*pi))*(self0less.*G1ygreat - self0great.*G1yless);
    Is0zbare = -(1/(2*pi))*(self0less.*G1zgreat - self0great.*G1zless);
    Is1xbare = 1i*(1/(2*pi))*(self1less.*G1xgreat - self1great.*G1xless);
    Is1ybare = -1i*(1/(2*pi))*(self1less.*G1ygreat - self1great.*G1yless);
    Is1zbare = -(1/(2*pi))*(self1less.*G1zgreat - self1great.*G1zless);
    Is0x = trapz(w,Is0xbare);
    Is0y = trapz(w,Is0ybare);
    Is0z = trapz(w,Is0zbare);
    Is1x = trapz(w,Is1xbare);
    Is1y = trapz(w,Is1ybare);
    Is1z = trapz(w,Is1zbare);
    Isx = Is0x + Is1x;
    Isy = Is0y + Is1y;
    Isz = Is0z + Is1z;

    %Energy current
    Ie0bare = -(1/(2*pi)).*w.*(self0less.*G0great - self0great.*G0less);
    Ie1bare = -(1/(2*pi)).*w.*(self1less.*G1zgreat - self1great.*G1zless);
    Ie0 = trapz(w,Ie0bare);
    Ie1 = trapz(w,Ie1bare);
    Ie = Ie0 + Ie1;

    %Heat current
    Iq = Ie - eV(1)*Ic;

    %In SI units
    %Iconv=1.602176565*10^(-19)/(6.58211928*10^(-16))*10^-3;%in A, elementary charge divided by hbar in eVs times 10^-3 as we use meV
    %Ic = Ic.*Iconv;
    %Isx = Isx.*Iconv;
    %Isy = Isy.*Iconv;
    %Isz = Isz.*Iconv;

    %energyconv = 1/(6.58211928*10^(-16))*10^-3;%1 divided by hbar in eVs times 10^-3 as we use meV
    %Ie = Ie*energyconv;
    %Iq = Iq*energyconv;

    DOS = 1i/(2*pi)*(G0great-G0less);
    MDOSz = 1i/(4*pi)*(G1zgreat-G1zless);
end
