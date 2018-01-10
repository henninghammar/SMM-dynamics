function [dSx, dSy, dSz] = spinequationofmotion(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, SxT, SyT, SzT, tau)

Sx = S(1);
Sy = S(2);
Sz = S(3);

dSxbare=jH.*(Sy.*SzT-Sz.*SyT)+Sy.*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz.*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy.*(jDMx.*SyT-jDMy.*SxT)+Sz.*(jDMz.*SxT-jDMx.*SzT);
dSybare=jH.*(Sz.*SxT-Sx.*SzT)+Sz.*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx.*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz.*(jDMy.*SzT-jDMz.*SyT)+Sx.*(jDMx.*SyT-jDMy.*SxT);
dSzbare=jH.*(Sx.*SyT-Sy.*SxT)+Sx.*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy.*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx.*(jDMz.*SxT-jDMx.*SzT)+Sy.*(jDMy.*SzT-jDMz.*SyT);

dSx=-Beffz.*Sy+Beffy.*Sz+trapz(tau,dSxbare);
dSy=-Beffx.*Sz+Beffz.*Sx+trapz(tau,dSybare);
dSz=-Beffy.*Sx+Beffx.*Sy+trapz(tau,dSzbare);

end
