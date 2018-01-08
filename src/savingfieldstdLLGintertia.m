%Fields with spin
SBx(j)=-Beffz.*Sy(j)+Beffy.*Sz(j);
SBy(j)=-Beffx.*Sz(j)+Beffz.*Sx(j);
SBz(j)=-Beffy.*Sx(j)+Beffx.*Sy(j);

Sejx(j)=ejz.*Sy(j)-ejy.*Sz(j);
Sejy(j)=ejx.*Sz(j)-ejz.*Sx(j);
Sejz(j)=ejy.*Sx(j)-ejx.*Sy(j);

Smx(j)=-mz.*Sy(j)+my.*Sz(j);
Smy(j)=-mx.*Sz(j)+mz.*Sx(j);
Smz(j)=-my.*Sx(j)+mx.*Sy(j);

%V�xelverkan Sxj*S
SxjHx(j)=trapz(tau,jH.*(Sy(j).*Sz(j)-Sz(j).*Sy(j)));
SxjHy(j)=trapz(tau,jH.*(Sz(j).*Sx(j)-Sx(j).*Sz(j)));
SxjHz(j)=trapz(tau,jH.*(Sx(j).*Sy(j)-Sy(j).*Sx(j)));
SxjIx(j)=trapz(tau,Sy(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j))-Sz(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j)));
SxjIy(j)=trapz(tau,Sz(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j))-Sx(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j)));
SxjIz(j)=trapz(tau,Sx(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j))-Sy(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j)));
SxjDMx(j)=trapz(tau,-Sy(j).*(jDMx.*Sy(j)-jDMy.*Sx(j))+Sz(j).*(jDMz.*Sx(j)-jDMx.*Sz(j)));
SxjDMy(j)=trapz(tau,-Sz(j).*(jDMy.*Sz(j)-jDMz.*Sy(j))+Sx(j).*(jDMx.*Sy(j)-jDMy.*Sx(j)));
SxjDMz(j)=trapz(tau,-Sx(j).*(jDMz.*Sx(j)-jDMx.*Sz(j))+Sy(j).*(jDMy.*Sz(j)-jDMz.*Sy(j)));

%Det effektiva magnetf�ltet best�ende av ej + B + J*S
dSxbare1=jH.*(Sy(j).*Sz(j)-Sz(j).*Sy(j))+Sy(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j))-Sz(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j))-Sy(j).*(jDMx.*Sy(j)-jDMy.*Sx(j))+Sz(j).*(jDMz.*Sx(j)-jDMx.*Sz(j));
dSybare1=jH.*(Sz(j).*Sx(j)-Sx(j).*Sz(j))+Sz(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j))-Sx(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j))-Sz(j).*(jDMy.*Sz(j)-jDMz.*Sy(j))+Sx(j).*(jDMx.*Sy(j)-jDMy.*Sx(j));
dSzbare1=jH.*(Sx(j).*Sy(j)-Sy(j).*Sx(j))+Sx(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j))-Sy(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j))-Sx(j).*(jDMz.*Sx(j)-jDMx.*Sz(j))+Sy(j).*(jDMy.*Sz(j)-jDMz.*Sy(j));
SBeffx(j)=-Beffz.*Sy(j)+Beffy.*Sz(j)+trapz(tau,dSxbare1);
SBeffy(j)=-Beffx.*Sz(j)+Beffz.*Sx(j)+trapz(tau,dSybare1);
SBeffz(j)=-Beffy.*Sx(j)+Beffx.*Sy(j)+trapz(tau,dSzbare1);

%Dessa parametrar ber�knar ut Gilbert-d�mpningen i de olika
%riktningarna
dSxbarejH=-(t(j)-tau).*jH.*(Sy(j).*dSz-Sz(j).*dSy);
dSybarejH=-(t(j)-tau).*jH.*(Sz(j).*dSx-Sx(j).*dSz);
dSzbarejH=-(t(j)-tau).*jH.*(Sx(j).*dSy-Sy(j).*dSx);
SjHx(j)=trapz(tau,dSxbarejH);
SjHy(j)=trapz(tau,dSybarejH);
SjHz(j)=trapz(tau,dSzbarejH);

dSxbarejDMx=(t(j)-tau).*(Sy(j).*(jDMx.*dSy-jDMy.*dSx)+Sz(j).*(jDMz.*dSx-jDMx.*dSz));
dSybarejDMy=(t(j)-tau).*(Sz(j).*(jDMy.*dSz-jDMz.*dSy)+Sx(j).*(jDMx.*dSy-jDMy.*dSx));
dSzbarejDMz=(t(j)-tau).*(Sx(j).*(jDMz.*dSx-jDMx.*dSz)+Sy(j).*(jDMy.*dSz-jDMz.*dSy));
SjDMx(j)=trapz(tau,dSxbarejDMx);
SjDMy(j)=trapz(tau,dSybarejDMy);
SjDMz(j)=trapz(tau,dSzbarejDMz);

dSxbarejIx=-(t(j)-tau).*(Sy(j).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-Sz(j).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz));
dSybarejIy=-(t(j)-tau).*(Sz(j).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-Sx(j).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz));
dSzbarejIz=-(t(j)-tau).*(Sx(j).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-Sy(j).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz));
SjIx(j)=trapz(tau,dSxbarejIx);
SjIy(j)=trapz(tau,dSybarejIy);
SjIz(j)=trapz(tau,dSzbarejIz);

%Dessa parametrar ber�knar ut tr�ghetstensorn i de olika
%riktningarna
ddSxbarejH=0.5*(t(j)-tau).^2.*jH.*(Sy(j).*ddSz-Sz(j).*ddSy);
ddSybarejH=0.5*(t(j)-tau).^2.*jH.*(Sz(j).*ddSx-Sx(j).*ddSz);
ddSzbarejH=0.5*(t(j)-tau).^2.*jH.*(Sx(j).*ddSy-Sy(j).*ddSx);
dSjHx(j)=trapz(tau,ddSxbarejH);
dSjHy(j)=trapz(tau,ddSybarejH);
dSjHz(j)=trapz(tau,ddSzbarejH);

ddSxbarejDMx=-0.5*(t(j)-tau).^2.*(Sy(j).*(jDMx.*ddSy-jDMy.*ddSx)+Sz(j).*(jDMz.*ddSx-jDMx.*ddSz));
ddSybarejDMy=-0.5*(t(j)-tau).^2.*(Sz(j).*(jDMy.*ddSz-jDMz.*ddSy)+Sx(j).*(jDMx.*ddSy-jDMy.*ddSx));
ddSzbarejDMz=-0.5*(t(j)-tau).^2.*(Sx(j).*(jDMz.*ddSx-jDMx.*ddSz)+Sy(j).*(jDMy.*ddSz-jDMz.*ddSy));
dSjDMx(j)=trapz(tau,ddSxbarejDMx);
dSjDMy(j)=trapz(tau,ddSybarejDMy);
dSjDMz(j)=trapz(tau,ddSzbarejDMz);

ddSxbarejIx=0.5*(t(j)-tau).^2.*(Sy(j).*(jIzx.*ddSx+jIzy.*dSy+jIzz.*ddSz)-Sz(j).*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz));
ddSybarejIy=0.5*(t(j)-tau).^2.*(Sz(j).*(jIxx.*ddSx+jIxy.*dSy+jIxz.*ddSz)-Sx(j).*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz));
ddSzbarejIz=0.5*(t(j)-tau).^2.*(Sx(j).*(jIyx.*ddSx+jIyy.*dSy+jIyz.*ddSz)-Sy(j).*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz));
dSjIx(j)=trapz(tau,ddSxbarejIx);
dSjIy(j)=trapz(tau,ddSybarejIy);
dSjIz(j)=trapz(tau,ddSzbarejIz);

%The fields acting on the spin
ejvect(:,j)=[ejx, ejy, ejz];

mvect(:,j)=[mx, my, mz];

Beffvect(:,j)=[Beffx, Beffy, Beffz];

jHt(j)=trapz(tau,jH);
jDMvect(:,j)=[trapz(tau,jDMx), trapz(tau,jDMy), trapz(tau,jDMz)];
jIvect(:,j)=[trapz(tau,jIxx), trapz(tau,jIyy), trapz(tau,jIzz), trapz(tau,jIxy), trapz(tau,jIxz), trapz(tau,jIyz)]
