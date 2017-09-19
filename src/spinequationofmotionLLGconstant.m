function [dSxnew, dSynew, dSznew] = spinequationofmotionLLGconstant(Beffx, Beffy, Beffz, jH, jIxx, jIyy, jIzz, jIxy, jIxz, jIyz, jIyx, jIzx, jIzy, jDMx, jDMy, jDMz, GjH, GjIxx, GjIyy, GjIzz, GjIxy, GjIxz, GjIyx, GjIzx, GjIzy, GjIyz, GjDMx, GjDMy, GjDMz, S, dS)

dSx = dS(1);
dSy = dS(2);
dSz = dS(3);

dSxbare1=jH.*(S(2).*S(3)-S(3).*S(2))+S(2).*(jIzx.*S(1)+jIzy.*S(2)+jIzz.*S(3))-S(3).*(jIyx.*S(1)+jIyy.*S(2)+jIyz.*S(3))-S(2).*(jDMx.*S(2)-jDMy.*S(1))+S(3).*(jDMz.*S(1)-jDMx.*S(3));
dSybare1=jH.*(S(3).*S(1)-S(1).*S(3))+S(3).*(jIxx.*S(1)+jIxy.*S(2)+jIxz.*S(3))-S(1).*(jIzx.*S(1)+jIzy.*S(2)+jIzz.*S(3))-S(3).*(jDMy.*S(3)-jDMz.*S(2))+S(1).*(jDMx.*S(2)-jDMy.*S(1));
dSzbare1=jH.*(S(1).*S(2)-S(2).*S(1))+S(1).*(jIyx.*S(1)+jIyy.*S(2)+jIyz.*S(3))-S(2).*(jIxx.*S(1)+jIxy.*S(2)+jIxz.*S(3))-S(1).*(jDMz.*S(1)-jDMx.*S(3))+S(2).*(jDMy.*S(3)-jDMz.*S(2));
dSxbare2=GjH.*(S(2).*dSz-S(3).*dSy)+S(2).*(GjIzx.*dSx+GjIzy.*dSy+GjIzz.*dSz)-S(3).*(GjIyx.*dSx+GjIyy.*dSy+GjIyz.*dSz)-S(2).*(GjDMx.*dSy-GjDMy.*dSx)+S(3).*(GjDMz.*dSx-GjDMx.*dSz);
dSybare2=GjH.*(S(3).*dSx-S(1).*dSz)+S(3).*(GjIxx.*dSx+GjIxy.*dSy+GjIxz.*dSz)-S(1).*(GjIzx.*dSx+GjIzy.*dSy+GjIzz.*dSz)-S(3).*(GjDMy.*dSz-GjDMz.*dSy)+S(1).*(GjDMx.*dSy-GjDMy.*dSx);
dSzbare2=GjH.*(S(1).*dSy-S(2).*dSx)+S(1).*(GjIyx.*dSx+GjIyy.*dSy+GjIyz.*dSz)-S(2).*(GjIxx.*dSx+GjIxy.*dSy+GjIxz.*dSz)-S(1).*(GjDMz.*dSx-GjDMx.*dSz)+S(2).*(GjDMy.*dSz-GjDMz.*dSy);
dSxnew=-Beffz.*S(2)+Beffy.*S(3)+dSxbare1+dSxbare2;
dSynew=-Beffx.*S(3)+Beffz.*S(1)+dSybare1+dSybare2;
dSznew=-Beffy.*S(1)+Beffx.*S(2)+dSzbare1+dSzbare2;

end
