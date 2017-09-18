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
dSx1=-Beffz.*Sy(j)+Beffy.*Sz(j)+trapz(tau,dSxbare1-(t(j)-tau).*dSxbare2);
dSy1=-Beffx.*Sz(j)+Beffz.*Sx(j)+trapz(tau,dSybare1-(t(j)-tau).*dSybare2);
dSz1=-Beffy.*Sx(j)+Beffx.*Sy(j)+trapz(tau,dSzbare1-(t(j)-tau).*dSzbare2);
