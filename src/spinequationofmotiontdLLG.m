function [dSx, dSy, dSz] = spinequationofmotiontdLLG(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, dS, tau, t)

dSxbare1=jH.*(S(2).*S(3)-S(3).*S(2))+S(2).*(jIzx.*S(1)+jIzy.*S(2)+jIzz.*S(3))-S(3).*(jIyx.*S(1)+jIyy.*S(2)+jIyz.*S(3))-S(2).*(jDMx.*S(2)-jDMy.*S(1))+S(3).*(jDMz.*S(1)-jDMx.*S(3));
dSybare1=jH.*(S(3).*S(1)-S(1).*S(3))+S(3).*(jIxx.*S(1)+jIxy.*S(2)+jIxz.*S(3))-S(1).*(jIzx.*S(1)+jIzy.*S(2)+jIzz.*S(3))-S(3).*(jDMy.*S(3)-jDMz.*S(2))+S(1).*(jDMx.*S(2)-jDMy.*S(1));
dSzbare1=jH.*(S(1).*S(2)-S(2).*S(1))+S(1).*(jIyx.*S(1)+jIyy.*S(2)+jIyz.*S(3))-S(2).*(jIxx.*S(1)+jIxy.*S(2)+jIxz.*S(3))-S(1).*(jDMz.*S(1)-jDMx.*S(3))+S(2).*(jDMy.*S(3)-jDMz.*S(2));
dSxbare2=jH.*(S(2).*dS(3)-S(3).*dS(2))+S(2).*(jIzx.*dS(1)+jIzy.*dS(2)+jIzz.*dS(3))-S(3).*(jIyx.*dS(1)+jIyy.*dS(2)+jIyz.*dS(3))-S(2).*(jDMx.*dS(2)-jDMy.*dS(1))+S(3).*(jDMz.*dS(1)-jDMx.*dS(3));
dSybare2=jH.*(S(3).*dS(1)-S(1).*dS(3))+S(3).*(jIxx.*dS(1)+jIxy.*dS(2)+jIxz.*dS(3))-S(1).*(jIzx.*dS(1)+jIzy.*dS(2)+jIzz.*dS(3))-S(3).*(jDMy.*dS(3)-jDMz.*dS(2))+S(1).*(jDMx.*dS(2)-jDMy.*dS(1));
dSzbare2=jH.*(S(1).*dS(2)-S(2).*dS(1))+S(1).*(jIyx.*dS(1)+jIyy.*dS(2)+jIyz.*dS(3))-S(2).*(jIxx.*dS(1)+jIxy.*dS(2)+jIxz.*dS(3))-S(1).*(jDMz.*dS(1)-jDMx.*dS(3))+S(2).*(jDMy.*dS(3)-jDMz.*dS(2));

dSx=-Beffz.*S(2)+Beffy.*S(3)+trapz(tau,dSxbare1-(t-tau).*dSxbare2);
dSy=-Beffx.*S(3)+Beffz.*S(1)+trapz(tau,dSybare1-(t-tau).*dSybare2);
dSz=-Beffy.*S(1)+Beffx.*S(2)+trapz(tau,dSzbare1-(t-tau).*dSzbare2);

end
