SBx(j)=-Beffz.*Sy(j)+Beffy.*Sz(j);
SBy(j)=-Beffx.*Sz(j)+Beffz.*Sx(j);
SBz(j)=-Beffy.*Sx(j)+Beffx.*Sy(j);

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
dSxbarejH(:,j)=-(t(j)-tau).*jH.*(Sy(j).*dSz-Sz(j).*dSy);
dSybarejH(:,j)=-(t(j)-tau).*jH.*(Sz(j).*dSx-Sx(j).*dSz);
dSzbarejH(:,j)=-(t(j)-tau).*jH.*(Sx(j).*dSy-Sy(j).*dSx);
SjHx(j)=trapz(tau,dSxbarejH(:,j));
SjHy(j)=trapz(tau,dSybarejH(:,j));
SjHz(j)=trapz(tau,dSzbarejH(:,j));

dSxbarejDMx(:,j)=(t(j)-tau).*(Sy(j).*(jDMx.*dSy-jDMy.*dSx)+Sz(j).*(jDMz.*dSx-jDMx.*dSz));
dSybarejDMy(:,j)=(t(j)-tau).*(Sz(j).*(jDMy.*dSz-jDMz.*dSy)+Sx(j).*(jDMx.*dSy-jDMy.*dSx));
dSzbarejDMz(:,j)=(t(j)-tau).*(Sx(j).*(jDMz.*dSx-jDMx.*dSz)+Sy(j).*(jDMy.*dSz-jDMz.*dSy));
SjDMx(j)=trapz(tau,dSxbarejDMx(:,j));
SjDMy(j)=trapz(tau,dSybarejDMy(:,j));
SjDMz(j)=trapz(tau,dSzbarejDMz(:,j));

dSxbarejIx(:,j)=-(t(j)-tau).*(Sy(j).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-Sz(j).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz));
dSybarejIy(:,j)=-(t(j)-tau).*(Sz(j).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-Sx(j).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz));
dSzbarejIz(:,j)=-(t(j)-tau).*(Sx(j).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-Sy(j).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz));
SjIx(j)=trapz(tau,dSxbarejIx(:,j));
SjIy(j)=trapz(tau,dSybarejIy(:,j));
SjIz(j)=trapz(tau,dSzbarejIz(:,j));

%The fields acting on the spin
Beffxvect(j)=Beffx;
Beffyvect(j)=Beffy;
Beffzvect(j)=Beffz;

barejHfield(:,j)=jH;
jHfield(j)=trapz(tau,barejHfield(:,j));

barejIfieldx(:,j)=jIxx.*dSx+jIxy.*dSy+jIxz.*dSz;
barejIfieldy(:,j)=jIyx.*dSx+jIyy.*dSy+jIyz.*dSz;
barejIfieldz(:,j)=jIzx.*dSx+jIzy.*dSy+jIzz.*dSz;
jIfieldx(j)=trapz(tau,barejIfieldx(:,j));
jIfieldy(j)=trapz(tau,barejIfieldy(:,j));
jIfieldz(j)=trapz(tau,barejIfieldz(:,j));

barejDMfieldx(:,j)=jDMz.*dSy-jDMy.*dSz;
barejDMfieldy(:,j)=jDMx.*dSz-jDMz.*dSx;
barejDMfieldz(:,j)=jDMy.*dSx-jDMx.*dSy;
jDMfieldx(j)=trapz(tau,barejDMfieldx(:,j));
jDMfieldy(j)=trapz(tau,barejDMfieldy(:,j));
jDMfieldz(j)=trapz(tau,barejDMfieldz(:,j));

jHtot(:,j)=jH;
jDMxtot(:,j)=jDMx;
jDMytot(:,j)=jDMy;
jDMztot(:,j)=jDMz;
jIxxtot(:,j)=jIxx;
jIyytot(:,j)=jIyy;
jIzztot(:,j)=jIzz;
jIxytot(:,j)=jIxy;
jIxztot(:,j)=jIxz;
jIyztot(:,j)=jIyz;

jHt(j)=trapz(tau,jH);
jDMxt(j)=trapz(tau,jDMx);
jDMyt(j)=trapz(tau,jDMy);
jDMzt(j)=trapz(tau,jDMz);
jIxxt(j)=trapz(tau,jIxx);
jIyyt(j)=trapz(tau,jIyy);
jIzzt(j)=trapz(tau,jIzz);
jIxyt(j)=trapz(tau,jIxy);
jIxzt(j)=trapz(tau,jIxz);
jIyzt(j)=trapz(tau,jIyz);
