%Changed such that it calculates the SEOM with Sx2, Sy2 och Sz2 och SxT2, SyT2, SzT2
SxT2=[SxT(2:end),Sx2];
SyT2=[SyT(2:end),Sy2];
SzT2=[SzT(2:end),Sz2];

Beffx = -ejx + J*mx;
Beffy = -ejy + J*my;
Beffz = wL - ejz + J*mz;

dSxbare=jH.*(Sy2.*SzT2-Sz2.*SyT2)+Sy2.*(jIzx.*SxT2+jIzy.*SyT2+jIzz.*SzT2)-Sz2.*(jIyx.*SxT2+jIyy.*SyT2+jIyz.*SzT2)-Sy2.*(jDMx.*SyT2-jDMy.*SxT2)+Sz2.*(jDMz.*SxT2-jDMx.*SzT2);
dSybare=jH.*(Sz2.*SxT2-Sx2.*SzT2)+Sz2.*(jIxx.*SxT2+jIxy.*SyT2+jIxz.*SzT2)-Sx2.*(jIzx.*SxT2+jIzy.*SyT2+jIzz.*SzT2)-Sz2.*(jDMy.*SzT2-jDMz.*SyT2)+Sx2.*(jDMx.*SyT2-jDMy.*SxT2);
dSzbare=jH.*(Sx2.*SyT2-Sy2.*SxT2)+Sx2.*(jIyx.*SxT2+jIyy.*SyT2+jIyz.*SzT2)-Sy2.*(jIxx.*SxT2+jIxy.*SyT2+jIxz.*SzT2)-Sx2.*(jDMz.*SxT2-jDMx.*SzT2)+Sy2.*(jDMy.*SzT2-jDMz.*SyT2);
dSx2=-Beffz.*Sy2+Beffy.*Sz2+trapz(tau,dSxbare);
dSy2=-Beffx.*Sz2+Beffz.*Sx2+trapz(tau,dSybare);
dSz2=-Beffy.*Sx2+Beffx.*Sy2+trapz(tau,dSzbare);
