%Calculation of particle and energy currents in the QD

InL0bare=KLgreat.*G0less+KLless.*G0great;
InL1bare=KLgreat.*G1zless+KLless.*G1zgreat;
InL0(j)=-4.*g0(1)*imag(trapz(tau,InL0bare));
InL1(j)=-4.*gS(1)*imag(trapz(tau,InL1bare));

InR0bare=KRgreat.*G0less+KRless.*G0great;
InR1bare=KRgreat.*G1zless+KRless.*G1zgreat;
InR0(j)=-4.*g0(2)*imag(trapz(tau,InR0bare));
InR1(j)=-4.*gS(2)*imag(trapz(tau,InR1bare));

IeL0bare=energyKLgreat.*G0less+energyKLless.*G0great;
IeL1bare=energyKLgreat.*G1zless+energyKLless.*G1zgreat;
IeL0(j)=-4.*g0(1)*imag(trapz(tau,IeL0bare));
IeL1(j)=-4.*gS(1)*imag(trapz(tau,IeL1bare));

IeR0bare=energyKRgreat.*G0less+energyKRless.*G0great;
IeR1bare=energyKRgreat.*G1zless+energyKRless.*G1zgreat;
IeR0(j)=-4.*g0(2)*imag(trapz(tau,IeR0bare));
IeR1(j)=-4.*gS(2)*imag(trapz(tau,IeR1bare));

IqL0(j)=IeL0(j)-mu(1)*InL0(j);
IqL1(j)=IeL1(j)-mu(1)*InL1(j);

IqR0(j)=IeR0(j)-mu(2)*InR0(j);
IqR1(j)=IeR1(j)-mu(2)*InR1(j);

% [G0lessent, G0greatent, G1xlessent, G1xgreatent, G1ylessent, G1ygreatent, G1zlessent, G1zgreatent] = greensfunction(t(j), t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);
% entropy0(j) = -imag(G0lessent).*log(imag(G0lessent));
% if G1zlessent == 0
%     entropy1(j) = 0;
% else 
%     entropy1(j) = -imag(G1zlessent).*log(imag(G1zlessent));
% end
% entropy(j) = entropy0(j) + entropy1(j);

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