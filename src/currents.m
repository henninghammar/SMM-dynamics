%Calculation of currents in the QD

Ic0bare=KLgreat.*G0less+KLless.*G0great;
Ic1bare=KLgreat.*G1zless+KLless.*G1zgreat;
Ic0(j)=-4.*g0(1)*imag(trapz(tau,Ic0bare));
Ic1(j)=-4.*gS(1)*imag(trapz(tau,Ic1bare));

Isx0bare=KLgreat.*G1xless+KLless.*G1xgreat;
Isx1bare=-1i.*(KLgreat.*G1yless+KLless.*G1ygreat);
Isx0(j)=-4.*g0(1)*imag(trapz(tau,Isx0bare));
Isx1(j)=-4.*gS(1)*imag(trapz(tau,Isx1bare));

Isy0bare=KLgreat.*G1yless+KLless.*G1ygreat;
Isy1bare=1i.*(KLgreat.*G1xless+KLless.*G1xgreat);
Isy0(j)=-4.*g0(1)*imag(trapz(tau,Isy0bare));
Isy1(j)=-4.*gS(1)*imag(trapz(tau,Isy1bare));

Isz0bare=KLgreat.*G1zless+KLless.*G1zgreat;
Isz1bare=KLgreat.*G0less+KLless.*G0great;
Isz0(j)=-4.*g0(1)*imag(trapz(tau,Isz0bare));
Isz1(j)=-4.*gS(1)*imag(trapz(tau,Isz1bare));
