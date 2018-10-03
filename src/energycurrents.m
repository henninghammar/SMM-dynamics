%Calculation of particle and energy currents in the QD

InL0bare=KLgreat.*G0less+KLless.*G0great;
InL1bare=KLgreat.*G1zless+KLless.*G1zgreat;
InL0(j)=4.*g0(1)*imag(trapz(tau,InL0bare));
InL1(j)=4.*gS(1)*imag(trapz(tau,InL1bare));

InR0bare=KRgreat.*G0less+KRless.*G0great;
InR1bare=KRgreat.*G1zless+KRless.*G1zgreat;
InR0(j)=4.*g0(2)*imag(trapz(tau,InR0bare));
InR1(j)=4.*gS(2)*imag(trapz(tau,InR1bare));

IeL0bare=energyKLgreat.*energyG0less+energyKLless.*energyG0great;
IeL1bare=energyKLgreat.*energyG1zless+energyKLless.*energyG1zgreat;
IeL0(j)=4.*g0(1)*imag(trapz(tauenergy,IeL0bare));
IeL1(j)=4.*gS(1)*imag(trapz(tauenergy,IeL1bare));

IeR0bare=energyKRgreat.*energyG0less+energyKRless.*energyG0great;
IeR1bare=energyKRgreat.*energyG1zless+energyKRless.*energyG1zgreat;
IeR0(j)=4.*g0(2)*imag(trapz(tauenergy,IeR0bare));
IeR1(j)=4.*gS(2)*imag(trapz(tauenergy,IeR1bare));

IqL0(j)=IeL0(j)-mu(1)*InL0(j);
IqL1(j)=IeL1(j)-mu(1)*InL1(j);

IqR0(j)=IeR0(j)-mu(2)*InR0(j);
IqR1(j)=IeR1(j)-mu(2)*InR1(j);

conv=1;%/(6.58211928*10^(-16))*10^-3;%in A, elementary charge divided by hbar in eVs times 10^-3
InL(j)=(InL0(j)+InL1(j)).*conv;
IeL(j)=(IeL0(j)+IeL1(j)).*conv;
IqL(j)=(IqL0(j)+IqL1(j)).*conv;
InR(j)=(InR0(j)+InR1(j)).*conv;
IeR(j)=(IeR0(j)+IeR1(j)).*conv;
IqR(j)=(IqR0(j)+IqR1(j)).*conv;

deltaIe(j) = IeL(j) + IeR(j);
deltaIq(j) = IqL(j) + IqR(j);
if j > 1
  thistory = [t0:tstep:t(j)];
  energytransfer(j) = trapz(thistory, deltaIe);
  heattransfer(j) = trapz(thistory, deltaIq);
else
  energytransfer(j) = deltaIe;
  heattransfer(j) = deltaIq;
end
