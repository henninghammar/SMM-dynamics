%Saving the interesting fields

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

dSxbarejH=jH.*(Sy(j).*SzT-Sz(j).*SyT);
dSybarejH=jH.*(Sz(j).*SxT-Sx(j).*SzT);
dSzbarejH=jH.*(Sx(j).*SyT-Sy(j).*SxT);
SjHx(j)=trapz(tau,dSxbarejH);
SjHy(j)=trapz(tau,dSybarejH);
SjHz(j)=trapz(tau,dSzbarejH);

dSxbarejDMx=-Sy(j).*(jDMx.*SyT-jDMy.*SxT)+Sz(j).*(jDMz.*SxT-jDMx.*SzT);
dSybarejDMy=-Sz(j).*(jDMy.*SzT-jDMz.*SyT)+Sx(j).*(jDMx.*SyT-jDMy.*SxT);
dSzbarejDMz=-Sx(j).*(jDMz.*SxT-jDMx.*SzT)+Sy(j).*(jDMy.*SzT-jDMz.*SyT);
SjDMx(j)=trapz(tau,dSxbarejDMx);
SjDMy(j)=trapz(tau,dSybarejDMy);
SjDMz(j)=trapz(tau,dSzbarejDMz);

dSxbarejIx=Sy(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT);
dSybarejIy=Sz(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT);
dSzbarejIz=Sx(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT);
SjIx(j)=trapz(tau,dSxbarejIx);
SjIy(j)=trapz(tau,dSybarejIy);
SjIz(j)=trapz(tau,dSzbarejIz);

%The fields acting on the spin
ejvect(:,j)=[ejx, ejy, ejz];

mvect(:,j)=[mx, my, mz];

Beffvect(:,j)=[Beffx, Beffy, Beffz];

jHt(j)=trapz(tau,jH);
jDMvect(:,j)=[trapz(tau,jDMx), trapz(tau,jDMy), trapz(tau,jDMz)];
jIvect(:,j)=[trapz(tau,jIxx), trapz(tau,jIyy), trapz(tau,jIzz), trapz(tau,jIxy), trapz(tau,jIxz), trapz(tau,jIyz)]
