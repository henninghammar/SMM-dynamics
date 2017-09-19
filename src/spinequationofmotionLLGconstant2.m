%The spin equation of motion for LLG solution
Beffx = -ejx + J*mx;
Beffy = -ejy + J*my;
Beffz = wL - ejz + J*mz;

dSxbare1=jH.*(Sy2.*Sz2-Sz2.*Sy2)+Sy2.*(jIzx.*Sx2+jIzy.*Sy2+jIzz.*Sz2)-Sz2.*(jIyx.*Sx2+jIyy.*Sy2+jIyz.*Sz2)-Sy2.*(jDMx.*Sy2-jDMy.*Sx2)+Sz2.*(jDMz.*Sx2-jDMx.*Sz2);
dSybare1=jH.*(Sz2.*Sx2-Sx2.*Sz2)+Sz2.*(jIxx.*Sx2+jIxy.*Sy2+jIxz.*Sz2)-Sx2.*(jIzx.*Sx2+jIzy.*Sy2+jIzz.*Sz2)-Sz2.*(jDMy.*Sz2-jDMz.*Sy2)+Sx2.*(jDMx.*Sy2-jDMy.*Sx2);
dSzbare1=jH.*(Sx2.*Sy2-Sy2.*Sx2)+Sx2.*(jIyx.*Sx2+jIyy.*Sy2+jIyz.*Sz2)-Sy2.*(jIxx.*Sx2+jIxy.*Sy2+jIxz.*Sz2)-Sx2.*(jDMz.*Sx2-jDMx.*Sz2)+Sy2.*(jDMy.*Sz2-jDMz.*Sy2);
dSxbare2=GjH.*(Sy2.*dSz-Sz2.*dSy)+Sy2.*(GjIzx.*dSx+GjIzy.*dSy+GjIzz.*dSz)-Sz2.*(GjIyx.*dSx+GjIyy.*dSy+GjIyz.*dSz)-Sy2.*(GjDMx.*dSy-GjDMy.*dSx)+Sz2.*(GjDMz.*dSx-GjDMx.*dSz);
dSybare2=GjH.*(Sz2.*dSx-Sx2.*dSz)+Sz2.*(GjIxx.*dSx+GjIxy.*dSy+GjIxz.*dSz)-Sx2.*(GjIzx.*dSx+GjIzy.*dSy+GjIzz.*dSz)-Sz2.*(GjDMy.*dSz-GjDMz.*dSy)+Sx2.*(GjDMx.*dSy-GjDMy.*dSx);
dSzbare2=GjH.*(Sx2.*dSy-Sy2.*dSx)+Sx2.*(GjIyx.*dSx+GjIyy.*dSy+GjIyz.*dSz)-Sy2.*(GjIxx.*dSx+GjIxy.*dSy+GjIxz.*dSz)-Sx2.*(GjDMz.*dSx-GjDMx.*dSz)+Sy2.*(GjDMy.*dSz-GjDMz.*dSy);
dSx2=-Beffz.*Sy2+Beffy.*Sz2+dSxbare1+dSxbare2;
dSy2=-Beffx.*Sz2+Beffz.*Sx2+dSybare1+dSybare2;
dSz2=-Beffy.*Sx2+Beffx.*Sy2+dSzbare1+dSzbare2;

dSx=(dSx1+dSx2)/2;
dSy=(dSy1+dSy2)/2;
dSz=(dSz1+dSz2)/2;
