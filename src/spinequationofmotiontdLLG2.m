%The spin equation of motion for LLG solution
Beffx = -ejx + J*mx;
Beffy = -ejy + J*my;
Beffz = wL - ejz + J*mz;

dSxbare1=jH.*(Sy2.*Sz2-Sz2.*Sy2)+Sy2.*(jIzx.*Sx2+jIzy.*Sy2+jIzz.*Sz2)-Sz2.*(jIyx.*Sx2+jIyy.*Sy2+jIyz.*Sz2)-Sy2.*(jDMx.*Sy2-jDMy.*Sx2)+Sz2.*(jDMz.*Sx2-jDMx.*Sz2);
dSybare1=jH.*(Sz2.*Sx2-Sx2.*Sz2)+Sz2.*(jIxx.*Sx2+jIxy.*Sy2+jIxz.*Sz2)-Sx2.*(jIzx.*Sx2+jIzy.*Sy2+jIzz.*Sz2)-Sz2.*(jDMy.*Sz2-jDMz.*Sy2)+Sx2.*(jDMx.*Sy2-jDMy.*Sx2);
dSzbare1=jH.*(Sx2.*Sy2-Sy2.*Sx2)+Sx2.*(jIyx.*Sx2+jIyy.*Sy2+jIyz.*Sz2)-Sy2.*(jIxx.*Sx2+jIxy.*Sy2+jIxz.*Sz2)-Sx2.*(jDMz.*Sx2-jDMx.*Sz2)+Sy2.*(jDMy.*Sz2-jDMz.*Sy2);
dSxbare2=jH.*(Sy2.*dSz1-Sz2.*dSy1)+Sy2.*(jIzx.*dSx1+jIzy.*dSy1+jIzz.*dSz1)-Sz2.*(jIyx.*dSx1+jIyy.*dSy1+jIyz.*dSz1)-Sy2.*(jDMx.*dSy1-jDMy.*dSx1)+Sz2.*(jDMz.*dSx1-jDMx.*dSz1);
dSybare2=jH.*(Sz2.*dSx1-Sx2.*dSz1)+Sz2.*(jIxx.*dSx1+jIxy.*dSy1+jIxz.*dSz1)-Sx2.*(jIzx.*dSx1+jIzy.*dSy1+jIzz.*dSz1)-Sz2.*(jDMy.*dSz1-jDMz.*dSy1)+Sx2.*(jDMx.*dSy1-jDMy.*dSx1);
dSzbare2=jH.*(Sx2.*dSy1-Sy2.*dSx1)+Sx2.*(jIyx.*dSx1+jIyy.*dSy1+jIyz.*dSz1)-Sy2.*(jIxx.*dSx1+jIxy.*dSy1+jIxz.*dSz1)-Sx2.*(jDMz.*dSx1-jDMx.*dSz1)+Sy2.*(jDMy.*dSz1-jDMz.*dSy1);
dSx2=-Beffz.*Sy2+Beffy.*Sz2+trapz(tau,dSxbare1-(t(j)-tau).*dSxbare2);
dSy2=-Beffx.*Sz2+Beffz.*Sx2+trapz(tau,dSybare1-(t(j)-tau).*dSybare2);
dSz2=-Beffy.*Sx2+Beffx.*Sy2+trapz(tau,dSzbare1-(t(j)-tau).*dSzbare2);

dSx=(dSx1+dSx2)/2;
dSy=(dSy1+dSy2)/2;
dSz=(dSz1+dSz2)/2;
