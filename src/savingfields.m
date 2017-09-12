%Saving the interesting fields 

%Fields with spin
SBx(j)=-Beffz.*Sy(j)+Beffy.*Sz(j);
SBy(j)=-Beffx.*Sz(j)+Beffz.*Sx(j);
SBz(j)=-Beffy.*Sx(j)+Beffx.*Sy(j);

dSxbarejH(:,j)=jH.*(Sy(j).*SzT-Sz(j).*SyT);
dSybarejH(:,j)=jH.*(Sz(j).*SxT-Sx(j).*SzT);
dSzbarejH(:,j)=jH.*(Sx(j).*SyT-Sy(j).*SxT);
SjHx(j)=trapz(tau,dSxbarejH(:,j));
SjHy(j)=trapz(tau,dSybarejH(:,j));
SjHz(j)=trapz(tau,dSzbarejH(:,j));

dSxbarejDMx(:,j)=-Sy(j).*(jDMx.*SyT-jDMy.*SxT)+Sz(j).*(jDMz.*SxT-jDMx.*SzT);
dSybarejDMy(:,j)=-Sz(j).*(jDMy.*SzT-jDMz.*SyT)+Sx(j).*(jDMx.*SyT-jDMy.*SxT);
dSzbarejDMz(:,j)=-Sx(j).*(jDMz.*SxT-jDMx.*SzT)+Sy(j).*(jDMy.*SzT-jDMz.*SyT);
SjDMx(j)=trapz(tau,dSxbarejDMx(:,j));
SjDMy(j)=trapz(tau,dSybarejDMy(:,j));
SjDMz(j)=trapz(tau,dSzbarejDMz(:,j));

dSxbarejIx(:,j)=Sy(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT);
dSybarejIy(:,j)=Sz(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT);
dSzbarejIz(:,j)=Sx(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT);
SjIx(j)=trapz(tau,dSxbarejIx(:,j));
SjIy(j)=trapz(tau,dSybarejIy(:,j));
SjIz(j)=trapz(tau,dSzbarejIz(:,j));     

%The fields acting on the spin
Beffxvect(j)=Beffx;
Beffyvect(j)=Beffy;
Beffzvect(j)=Beffz;

barejHfield(:,j)=jH;
jHfield(j)=trapz(tau,barejHfield(:,j));

barejIfieldx(:,j)=jIxx.*SxT+jIxy.*SyT+jIxz.*SzT;
barejIfieldy(:,j)=jIyx.*SxT+jIyy.*SyT+jIyz.*SzT;
barejIfieldz(:,j)=jIzx.*SxT+jIzy.*SyT+jIzz.*SzT;
jIfieldx(j)=trapz(tau,barejIfieldx(:,j));
jIfieldy(j)=trapz(tau,barejIfieldy(:,j));
jIfieldz(j)=trapz(tau,barejIfieldz(:,j));  

barejDMfieldx(:,j)=jDMz.*SyT-jDMy.*SzT;
barejDMfieldy(:,j)=jDMx.*SzT-jDMz.*SxT;
barejDMfieldz(:,j)=jDMy.*SxT-jDMx.*SyT;
jDMfieldx(j)=trapz(tau,barejDMfieldx(:,j));
jDMfieldy(j)=trapz(tau,barejDMfieldy(:,j));
jDMfieldz(j)=trapz(tau,barejDMfieldz(:,j)); 

jHtot(:,j)=jH;
jDMxtot(:,j)=jDMx;
jDMytot(:,j)=jDMy;    
jDMztot(:,j)=jDMz;
jIxxtot(:,j)=jIxx;
jIyytot(:,j)=jIyy;
jIzztot(:,j)=jIzz;
jIxytot(:,j)=jIxy;
jIxztot(:,j)=jIxz;
jIyztot(:,j)=jIyz;

jHt(j)=trapz(tau,jH);
jDMxt(j)=trapz(tau,jDMx);
jDMyt(j)=trapz(tau,jDMy);    
jDMzt(j)=trapz(tau,jDMz);
jIxxt(j)=trapz(tau,jIxx);
jIyyt(j)=trapz(tau,jIyy);
jIzzt(j)=trapz(tau,jIzz);
jIxyt(j)=trapz(tau,jIxy);
jIxzt(j)=trapz(tau,jIxz);
jIyzt(j)=trapz(tau,jIyz);