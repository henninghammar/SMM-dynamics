function [dSx, dSy, dSz] = spinequationofmotion(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, ST, tau)

dSxbare=jH.*(S(2).*ST(3)-S(3).*ST(2))+S(2).*(jIzx.*ST(1)+jIzy.*ST(2)+jIzz.*ST(3))-S(3).*(jIyx.*ST(1)+jIyy.*ST(2)+jIyz.*ST(3))-S(2).*(jDMx.*ST(2)-jDMy.*ST(1))+S(3).*(jDMz.*ST(1)-jDMx.*ST(3));
dSybare=jH.*(S(3).*ST(1)-S(1).*ST(3))+S(3).*(jIxx.*ST(1)+jIxy.*ST(2)+jIxz.*ST(3))-S(1).*(jIzx.*ST(1)+jIzy.*ST(2)+jIzz.*ST(3))-S(3).*(jDMy.*ST(3)-jDMz.*ST(2))+S(1).*(jDMx.*ST(2)-jDMy.*ST(1));
dSzbare=jH.*(S(1).*ST(2)-S(2).*ST(1))+S(1).*(jIyx.*ST(1)+jIyy.*ST(2)+jIyz.*ST(3))-S(2).*(jIxx.*ST(1)+jIxy.*ST(2)+jIxz.*ST(3))-S(1).*(jDMz.*ST(1)-jDMx.*ST(3))+S(2).*(jDMy.*ST(3)-jDMz.*ST(2));

dSx=-Beffz.*S(2)+Beffy.*S(3)+trapz(tau,dSxbare);
dSy=-Beffx.*S(3)+Beffz.*S(1)+trapz(tau,dSybare);
dSz=-Beffy.*S(1)+Beffx.*S(2)+trapz(tau,dSzbare);

end
