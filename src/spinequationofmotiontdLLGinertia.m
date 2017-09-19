function [dSxnew, dSynew, dSznew, ddSxnew, ddSynew, ddSznew] = spinequationofmotiontdLLGinertia(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, dS, ddS, tau, t)

dSx = dS(1);
dSy = dS(2);
dSz = dS(3);
ddSx = ddS(1);
ddSy = ddS(2);
ddSz = ddS(3);

dSxbare1=jH.*(S(2).*S(3)-S(3).*S(2))+S(2).*(jIzx.*S(1)+jIzy.*S(2)+jIzz.*S(3))-S(3).*(jIyx.*S(1)+jIyy.*S(2)+jIyz.*S(3))-S(2).*(jDMx.*S(2)-jDMy.*S(1))+S(3).*(jDMz.*S(1)-jDMx.*S(3));
dSybare1=jH.*(S(3).*S(1)-S(1).*S(3))+S(3).*(jIxx.*S(1)+jIxy.*S(2)+jIxz.*S(3))-S(1).*(jIzx.*S(1)+jIzy.*S(2)+jIzz.*S(3))-S(3).*(jDMy.*S(3)-jDMz.*S(2))+S(1).*(jDMx.*S(2)-jDMy.*S(1));
dSzbare1=jH.*(S(1).*S(2)-S(2).*S(1))+S(1).*(jIyx.*S(1)+jIyy.*S(2)+jIyz.*S(3))-S(2).*(jIxx.*S(1)+jIxy.*S(2)+jIxz.*S(3))-S(1).*(jDMz.*S(1)-jDMx.*S(3))+S(2).*(jDMy.*S(3)-jDMz.*S(2));
dSxbare2=jH.*(S(2).*dSz-S(3).*dSy)+S(2).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-S(3).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-S(2).*(jDMx.*dSy-jDMy.*dSx)+S(3).*(jDMz.*dSx-jDMx.*dSz);
dSybare2=jH.*(S(3).*dSx-S(1).*dSz)+S(3).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-S(1).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-S(3).*(jDMy.*dSz-jDMz.*dSy)+S(1).*(jDMx.*dSy-jDMy.*dSx);
dSzbare2=jH.*(S(1).*dSy-S(2).*dSx)+S(1).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-S(2).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-S(1).*(jDMz.*dSx-jDMx.*dSz)+S(2).*(jDMy.*dSz-jDMz.*dSy);
dSxbare3=jH.*(S(2).*ddSz-S(3).*ddSy)+S(2).*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-S(3).*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-S(2).*(jDMx.*ddSy-jDMy.*ddSx)+S(3).*(jDMz.*ddSx-jDMx.*ddSz);
dSybare3=jH.*(S(3).*ddSx-S(1).*ddSz)+S(3).*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-S(1).*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-S(3).*(jDMy.*ddSz-jDMz.*ddSy)+S(1).*(jDMx.*ddSy-jDMy.*ddSx);
dSzbare3=jH.*(S(1).*ddSy-S(2).*ddSx)+S(1).*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-S(2).*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-S(1).*(jDMz.*ddSx-jDMx.*ddSz)+S(2).*(jDMy.*ddSz-jDMz.*ddSy);
dSxbare4=jH.*(dSy.*dSz-dSz.*dSy)+dSy.*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-dSz.*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-dSy.*(jDMx.*dSy-jDMy.*dSx)+dSz.*(jDMz.*dSx-jDMx.*dSz);
dSybare4=jH.*(dSz.*dSx-dSx.*dSz)+dSz.*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-dSx.*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-dSz.*(jDMy.*dSz-jDMz.*dSy)+dSx.*(jDMx.*dSy-jDMy.*dSx);
dSzbare4=jH.*(dSx.*dSy-dSy.*dSx)+dSx.*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-dSy.*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-dSx.*(jDMz.*dSx-jDMx.*dSz)+dSy.*(jDMy.*dSz-jDMz.*dSy);
dSxbare5=jH.*(dSy.*ddSz-dSz.*ddSy)+dSy.*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-dSz.*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-dSy.*(jDMx.*ddSy-jDMy.*ddSx)+dSz.*(jDMz.*ddSx-jDMx.*ddSz);
dSybare5=jH.*(dSz.*ddSx-dSx.*ddSz)+dSz.*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-dSx.*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-dSz.*(jDMy.*ddSz-jDMz.*ddSy)+dSx.*(jDMx.*ddSy-jDMy.*ddSx);
dSzbare5=jH.*(dSx.*ddSy-dSy.*ddSx)+dSx.*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-dSy.*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-dSx.*(jDMz.*ddSx-jDMx.*ddSz)+dSy.*(jDMy.*ddSz-jDMz.*ddSy);
dSxnew=-Beffz.*S(2)+Beffy.*S(3)+trapz(tau,dSxbare1-(t-tau).*dSxbare2+0.5*(t-tau).^2.*dSxbare3);
dSynew=-Beffx.*S(3)+Beffz.*S(1)+trapz(tau,dSybare1-(t-tau).*dSybare2+0.5*(t-tau).^2.*dSybare3);
dSznew=-Beffy.*S(1)+Beffx.*S(2)+trapz(tau,dSzbare1-(t-tau).*dSzbare2+0.5*(t-tau).^2.*dSzbare3);
ddSxnew=-Beffz.*S(2)+Beffy.*S(3)+trapz(tau,dSxbare4-(t-tau).*dSxbare5);
ddSynew=-Beffx.*S(3)+Beffz.*S(1)+trapz(tau,dSybare4-(t-tau).*dSybare5);
ddSznew=-Beffy.*S(1)+Beffx.*S(2)+trapz(tau,dSzbare4-(t-tau).*dSzbare5);

end
