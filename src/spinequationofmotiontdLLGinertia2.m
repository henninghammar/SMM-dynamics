%The spin equation of motion for LLG solution
SxT2=[SxT(2:end),Sx2];
SyT2=[SyT(2:end),Sy2];
SzT2=[SzT(2:end),Sz2];

Beffx = -ejx + J*mx;
Beffy = -ejy + J*my;
Beffz = wL - ejz + J*mz;

dSxbare1=jH.*(Sy2.*Sz2-Sz2.*Sy2)+Sy2.*(jIzx.*Sx2+jIzy.*Sy2+jIzz.*Sz2)-Sz2.*(jIyx.*Sx2+jIyy.*Sy2+jIyz.*Sz2)-Sy2.*(jDMx.*Sy2-jDMy.*Sx2)+Sz2.*(jDMz.*Sx2-jDMx.*Sz2);
dSybare1=jH.*(Sz2.*Sx2-Sx2.*Sz2)+Sz2.*(jIxx.*Sx2+jIxy.*Sy2+jIxz.*Sz2)-Sx2.*(jIzx.*Sx2+jIzy.*Sy2+jIzz.*Sz2)-Sz2.*(jDMy.*Sz2-jDMz.*Sy2)+Sx2.*(jDMx.*Sy2-jDMy.*Sx2);
dSzbare1=jH.*(Sx2.*Sy2-Sy2.*Sx2)+Sx2.*(jIyx.*Sx2+jIyy.*Sy2+jIyz.*Sz2)-Sy2.*(jIxx.*Sx2+jIxy.*Sy2+jIxz.*Sz2)-Sx2.*(jDMz.*Sx2-jDMx.*Sz2)+Sy2.*(jDMy.*Sz2-jDMz.*Sy2);
dSxbare2=jH.*(Sy2.*dSz1-Sz2.*dSy1)+Sy2.*(jIzx.*dSx1+jIzy.*dSy1+jIzz.*dSz1)-Sz2.*(jIyx.*dSx1+jIyy.*dSy1+jIyz.*dSz1)-Sy2.*(jDMx.*dSy1-jDMy.*dSx1)+Sz2.*(jDMz.*dSx1-jDMx.*dSz1);
dSybare2=jH.*(Sz2.*dSx1-Sx2.*dSz1)+Sz2.*(jIxx.*dSx1+jIxy.*dSy1+jIxz.*dSz1)-Sx2.*(jIzx.*dSx1+jIzy.*dSy1+jIzz.*dSz1)-Sz2.*(jDMy.*dSz1-jDMz.*dSy1)+Sx2.*(jDMx.*dSy1-jDMy.*dSx1);
dSzbare2=jH.*(Sx2.*dSy1-Sy2.*dSx1)+Sx2.*(jIyx.*dSx1+jIyy.*dSy1+jIyz.*dSz1)-Sy2.*(jIxx.*dSx1+jIxy.*dSy1+jIxz.*dSz1)-Sx2.*(jDMz.*dSx1-jDMx.*dSz1)+Sy2.*(jDMy.*dSz1-jDMz.*dSy1);
dSxbare3=jH.*(Sy2.*ddSz1-Sz2.*ddSy1)+Sy2.*(jIzx.*ddSx1+jIzy.*ddSy1+jIzz.*ddSz1)-Sz2.*(jIyx.*ddSx1+jIyy.*ddSy1+jIyz.*ddSz1)-Sy2.*(jDMx.*ddSy1-jDMy.*ddSx1)+Sz2.*(jDMz.*ddSx1-jDMx.*ddSz1);
dSybare3=jH.*(Sz2.*ddSx1-Sx2.*ddSz1)+Sz2.*(jIxx.*ddSx1+jIxy.*ddSy1+jIxz.*ddSz1)-Sx2.*(jIzx.*ddSx1+jIzy.*ddSy1+jIzz.*ddSz1)-Sz2.*(jDMy.*ddSz1-jDMz.*ddSy1)+Sx2.*(jDMx.*ddSy1-jDMy.*ddSx1);
dSzbare3=jH.*(Sx2.*ddSy1-Sy2.*ddSx1)+Sx2.*(jIyx.*ddSx1+jIyy.*ddSy1+jIyz.*ddSz1)-Sy2.*(jIxx.*ddSx1+jIxy.*ddSy1+jIxz.*ddSz1)-Sx2.*(jDMz.*ddSx1-jDMx.*ddSz1)+Sy2.*(jDMy.*ddSz1-jDMz.*ddSy1);
dSxbare4=jH.*(dSy1.*dSz1-dSz1.*dSy1)+dSy1.*(jIzx.*dSx1+jIzy.*dSy1+jIzz.*dSz1)-dSz1.*(jIyx.*dSx1+jIyy.*dSy1+jIyz.*dSz1)-dSy1.*(jDMx.*dSy1-jDMy.*dSx1)+dSz1.*(jDMz.*dSx1-jDMx.*dSz1);
dSybare4=jH.*(dSz1.*dSx1-dSx1.*dSz1)+dSz1.*(jIxx.*dSx1+jIxy.*dSy1+jIxz.*dSz1)-dSx1.*(jIzx.*dSx1+jIzy.*dSy1+jIzz.*dSz1)-dSz1.*(jDMy.*dSz1-jDMz.*dSy1)+dSx1.*(jDMx.*dSy1-jDMy.*dSx1);
dSzbare4=jH.*(dSx1.*dSy1-dSy1.*dSx1)+dSx1.*(jIyx.*dSx1+jIyy.*dSy1+jIyz.*dSz1)-dSy1.*(jIxx.*dSx1+jIxy.*dSy1+jIxz.*dSz1)-dSx1.*(jDMz.*dSx1-jDMx.*dSz1)+dSy1.*(jDMy.*dSz1-jDMz.*dSy1);
dSxbare5=jH.*(dSy1.*ddSz1-dSz1.*ddSy1)+dSy1.*(jIzx.*ddSx1+jIzy.*ddSy1+jIzz.*ddSz1)-dSz1.*(jIyx.*ddSx1+jIyy.*ddSy1+jIyz.*ddSz1)-dSy1.*(jDMx.*ddSy1-jDMy.*ddSx1)+dSz1.*(jDMz.*ddSx1-jDMx.*ddSz1);
dSybare5=jH.*(dSz1.*ddSx1-dSx1.*ddSz1)+dSz1.*(jIxx.*ddSx1+jIxy.*ddSy1+jIxz.*ddSz1)-dSx1.*(jIzx.*ddSx1+jIzy.*ddSy1+jIzz.*ddSz1)-dSz1.*(jDMy.*ddSz1-jDMz.*ddSy1)+dSx1.*(jDMx.*ddSy1-jDMy.*ddSx1);
dSzbare5=jH.*(dSx1.*ddSy1-dSy1.*ddSx1)+dSx1.*(jIyx.*ddSx1+jIyy.*ddSy1+jIyz.*ddSz1)-dSy1.*(jIxx.*ddSx1+jIxy.*ddSy1+jIxz.*ddSz1)-dSx1.*(jDMz.*ddSx1-jDMx.*ddSz1)+dSy1.*(jDMy.*ddSz1-jDMz.*ddSy1);
dSx2=-Beffz.*Sy2+Beffy.*Sz2+trapz(tau,dSxbare1-(t(j)-tau).*dSxbare2+0.5*(t(j)-tau).^2.*dSxbare3);
dSy2=-Beffx.*Sz2+Beffz.*Sx2+trapz(tau,dSybare1-(t(j)-tau).*dSybare2+0.5*(t(j)-tau).^2.*dSybare3);
dSz2=-Beffy.*Sx2+Beffx.*Sy2+trapz(tau,dSzbare1-(t(j)-tau).*dSzbare2+0.5*(t(j)-tau).^2.*dSzbare3);
ddSx2=-Beffz.*Sy2+Beffy.*Sz2+trapz(tau,dSxbare4-(t(j)-tau).*dSxbare5);
ddSy2=-Beffx.*Sz2+Beffz.*Sx2+trapz(tau,dSybare4-(t(j)-tau).*dSybare5);
ddSz2=-Beffy.*Sx2+Beffx.*Sy2+trapz(tau,dSzbare4-(t(j)-tau).*dSzbare5);

dSx=(dSx1+dSx2)/2;
dSy=(dSy1+dSy2)/2;
dSz=(dSz1+dSz2)/2;

ddSx=(ddSx1+ddSx2)/2;
ddSy=(ddSy1+ddSy2)/2;
ddSz=(ddSz1+ddSz2)/2;
