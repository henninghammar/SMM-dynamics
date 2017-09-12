%Calculating the spin equation of motion

Beffx = -ejx;
Beffy = -ejy;
Beffz = wL-ejz;

dSxbare=jH.*(Sy(j).*SzT-Sz(j).*SyT)+Sy(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy(j).*(jDMx.*SyT-jDMy.*SxT)+Sz(j).*(jDMz.*SxT-jDMx.*SzT);
dSybare=jH.*(Sz(j).*SxT-Sx(j).*SzT)+Sz(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz(j).*(jDMy.*SzT-jDMz.*SyT)+Sx(j).*(jDMx.*SyT-jDMy.*SxT);
dSzbare=jH.*(Sx(j).*SyT-Sy(j).*SxT)+Sx(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx(j).*(jDMz.*SxT-jDMx.*SzT)+Sy(j).*(jDMy.*SzT-jDMz.*SyT);
dSx1=-Beffz.*Sy(j)+Beffy.*Sz(j)+trapz(tau,dSxbare);
dSy1=-Beffx.*Sz(j)+Beffz.*Sx(j)+trapz(tau,dSybare);
dSz1=-Beffy.*Sx(j)+Beffx.*Sy(j)+trapz(tau,dSzbare);
