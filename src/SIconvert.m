%Converting time and currents to SI-units
tconv=6.5821189*10^(-16)/10^-3;%in s, hbar in eVs/meV
t=tconv.*t;
Iconv=1.602176565*10^(-19)/(6.58211928*10^(-16))*10^-3;%in A, elementary charge divided by hbar in eVs times 10^-3

Ic=(Ic0+Ic1).*Iconv;
Isx=(Isx0+Isx1).*Iconv;
Isy=(Isy0+Isy1).*Iconv;
Isz=(Isz0+Isz1).*Iconv;