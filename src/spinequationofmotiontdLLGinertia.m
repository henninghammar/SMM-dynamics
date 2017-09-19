%The spin equation of motion for LLG solution
if j == 1
   SxT=-Sxy*sin(wL*tau);
   SyT=Sxy*cos(wL*tau);
   SzT=Sz0*ones(1,length(tau));
else
   SxT=[SxT(2:end),Sx(j)];
   SyT=[SyT(2:end),Sy(j)];
   SzT=[SzT(2:end),Sz(j)];
end

Beffx = -ejx + J*mx;
Beffy = -ejy + J*my;
Beffz = wL - ejz + J*mz;

dSxbare1=jH.*(Sy(j).*Sz(j)-Sz(j).*Sy(j))+Sy(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j))-Sz(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j))-Sy(j).*(jDMx.*Sy(j)-jDMy.*Sx(j))+Sz(j).*(jDMz.*Sx(j)-jDMx.*Sz(j));
dSybare1=jH.*(Sz(j).*Sx(j)-Sx(j).*Sz(j))+Sz(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j))-Sx(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j))-Sz(j).*(jDMy.*Sz(j)-jDMz.*Sy(j))+Sx(j).*(jDMx.*Sy(j)-jDMy.*Sx(j));
dSzbare1=jH.*(Sx(j).*Sy(j)-Sy(j).*Sx(j))+Sx(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j))-Sy(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j))-Sx(j).*(jDMz.*Sx(j)-jDMx.*Sz(j))+Sy(j).*(jDMy.*Sz(j)-jDMz.*Sy(j));
dSxbare2=jH.*(Sy(j).*dSz-Sz(j).*dSy)+Sy(j).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-Sz(j).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-Sy(j).*(jDMx.*dSy-jDMy.*dSx)+Sz(j).*(jDMz.*dSx-jDMx.*dSz);
dSybare2=jH.*(Sz(j).*dSx-Sx(j).*dSz)+Sz(j).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-Sx(j).*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-Sz(j).*(jDMy.*dSz-jDMz.*dSy)+Sx(j).*(jDMx.*dSy-jDMy.*dSx);
dSzbare2=jH.*(Sx(j).*dSy-Sy(j).*dSx)+Sx(j).*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-Sy(j).*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-Sx(j).*(jDMz.*dSx-jDMx.*dSz)+Sy(j).*(jDMy.*dSz-jDMz.*dSy);
dSxbare3=jH.*(Sy(j).*ddSz-Sz(j).*ddSy)+Sy(j).*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-Sz(j).*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-Sy(j).*(jDMx.*ddSy-jDMy.*ddSx)+Sz(j).*(jDMz.*ddSx-jDMx.*ddSz);
dSybare3=jH.*(Sz(j).*ddSx-Sx(j).*ddSz)+Sz(j).*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-Sx(j).*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-Sz(j).*(jDMy.*ddSz-jDMz.*ddSy)+Sx(j).*(jDMx.*ddSy-jDMy.*ddSx);
dSzbare3=jH.*(Sx(j).*ddSy-Sy(j).*ddSx)+Sx(j).*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-Sy(j).*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-Sx(j).*(jDMz.*ddSx-jDMx.*ddSz)+Sy(j).*(jDMy.*ddSz-jDMz.*ddSy);
dSxbare4=jH.*(dSy.*dSz-dSz.*dSy)+dSy.*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-dSz.*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-dSy.*(jDMx.*dSy-jDMy.*dSx)+dSz.*(jDMz.*dSx-jDMx.*dSz);
dSybare4=jH.*(dSz.*dSx-dSx.*dSz)+dSz.*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-dSx.*(jIzx.*dSx+jIzy.*dSy+jIzz.*dSz)-dSz.*(jDMy.*dSz-jDMz.*dSy)+dSx.*(jDMx.*dSy-jDMy.*dSx);
dSzbare4=jH.*(dSx.*dSy-dSy.*dSx)+dSx.*(jIyx.*dSx+jIyy.*dSy+jIyz.*dSz)-dSy.*(jIxx.*dSx+jIxy.*dSy+jIxz.*dSz)-dSx.*(jDMz.*dSx-jDMx.*dSz)+dSy.*(jDMy.*dSz-jDMz.*dSy);
dSxbare5=jH.*(dSy.*ddSz-dSz.*ddSy)+dSy.*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-dSz.*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-dSy.*(jDMx.*ddSy-jDMy.*ddSx)+dSz.*(jDMz.*ddSx-jDMx.*ddSz);
dSybare5=jH.*(dSz.*ddSx-dSx.*ddSz)+dSz.*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-dSx.*(jIzx.*ddSx+jIzy.*ddSy+jIzz.*ddSz)-dSz.*(jDMy.*ddSz-jDMz.*ddSy)+dSx.*(jDMx.*ddSy-jDMy.*ddSx);
dSzbare5=jH.*(dSx.*ddSy-dSy.*ddSx)+dSx.*(jIyx.*ddSx+jIyy.*ddSy+jIyz.*ddSz)-dSy.*(jIxx.*ddSx+jIxy.*ddSy+jIxz.*ddSz)-dSx.*(jDMz.*ddSx-jDMx.*ddSz)+dSy.*(jDMy.*ddSz-jDMz.*ddSy);
dSx1=-Beffz.*Sy(j)+Beffy.*Sz(j)+trapz(tau,dSxbare1-(t(j)-tau).*dSxbare2+0.5*(t(j)-tau).^2.*dSxbare3);
dSy1=-Beffx.*Sz(j)+Beffz.*Sx(j)+trapz(tau,dSybare1-(t(j)-tau).*dSybare2+0.5*(t(j)-tau).^2.*dSybare3);
dSz1=-Beffy.*Sx(j)+Beffx.*Sy(j)+trapz(tau,dSzbare1-(t(j)-tau).*dSzbare2+0.5*(t(j)-tau).^2.*dSzbare3);
ddSx1=-Beffz.*Sy(j)+Beffy.*Sz(j)+trapz(tau,dSxbare4-(t(j)-tau).*dSxbare5);
ddSy1=-Beffx.*Sz(j)+Beffz.*Sx(j)+trapz(tau,dSybare4-(t(j)-tau).*dSybare5);
ddSz1=-Beffy.*Sx(j)+Beffx.*Sy(j)+trapz(tau,dSzbare4-(t(j)-tau).*dSzbare5);
