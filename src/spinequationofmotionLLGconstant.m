%The spin equation of motion for LLG solution
Beffx = -ejx + J*mx;
Beffy = -ejy + J*my;
Beffz = wL - ejz + J*mz;

dSxbare1=jH.*(Sy(j).*Sz(j)-Sz(j).*Sy(j))+Sy(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j))-Sz(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j))-Sy(j).*(jDMx.*Sy(j)-jDMy.*Sx(j))+Sz(j).*(jDMz.*Sx(j)-jDMx.*Sz(j));
dSybare1=jH.*(Sz(j).*Sx(j)-Sx(j).*Sz(j))+Sz(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j))-Sx(j).*(jIzx.*Sx(j)+jIzy.*Sy(j)+jIzz.*Sz(j))-Sz(j).*(jDMy.*Sz(j)-jDMz.*Sy(j))+Sx(j).*(jDMx.*Sy(j)-jDMy.*Sx(j));
dSzbare1=jH.*(Sx(j).*Sy(j)-Sy(j).*Sx(j))+Sx(j).*(jIyx.*Sx(j)+jIyy.*Sy(j)+jIyz.*Sz(j))-Sy(j).*(jIxx.*Sx(j)+jIxy.*Sy(j)+jIxz.*Sz(j))-Sx(j).*(jDMz.*Sx(j)-jDMx.*Sz(j))+Sy(j).*(jDMy.*Sz(j)-jDMz.*Sy(j));
dSxbare2=GjH.*(Sy(j).*dSz-Sz(j).*dSy)+Sy(j).*(GjIzx.*dSx+GjIzy.*dSy+GjIzz.*dSz)-Sz(j).*(GjIyx.*dSx+GjIyy.*dSy+GjIyz.*dSz)-Sy(j).*(GjDMx.*dSy-GjDMy.*dSx)+Sz(j).*(GjDMz.*dSx-GjDMx.*dSz);
dSybare2=GjH.*(Sz(j).*dSx-Sx(j).*dSz)+Sz(j).*(GjIxx.*dSx+GjIxy.*dSy+GjIxz.*dSz)-Sx(j).*(GjIzx.*dSx+GjIzy.*dSy+GjIzz.*dSz)-Sz(j).*(GjDMy.*dSz-GjDMz.*dSy)+Sx(j).*(GjDMx.*dSy-GjDMy.*dSx);
dSzbare2=GjH.*(Sx(j).*dSy-Sy(j).*dSx)+Sx(j).*(GjIyx.*dSx+GjIyy.*dSy+GjIyz.*dSz)-Sy(j).*(GjIxx.*dSx+GjIxy.*dSy+GjIxz.*dSz)-Sx(j).*(GjDMz.*dSx-GjDMx.*dSz)+Sy(j).*(GjDMy.*dSz-GjDMz.*dSy);
dSx1=-Beffz.*Sy(j)+Beffy.*Sz(j)+dSxbare1+dSxbare2;
dSy1=-Beffx.*Sz(j)+Beffz.*Sx(j)+dSybare1+dSybare2;
dSz1=-Beffy.*Sx(j)+Beffx.*Sy(j)+dSzbare1+dSzbare2;
