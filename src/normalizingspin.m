%Calculate new spin and normalizing it
dSx=(dSx1+dSx2)/2;
dSy=(dSy1+dSy2)/2;
dSz=(dSz1+dSz2)/2;

dSxvect1(j)=dSx1;
dSyvect1(j)=dSy1;
dSzvect1(j)=dSz1;
dSxvect2(j)=dSx2;
dSyvect2(j)=dSy2;
dSzvect2(j)=dSz2;
dSxvectTot(j)=dSx;
dSyvectTot(j)=dSy;
dSzvectTot(j)=dSz;

if j < length(t)
    Sx(j+1)=Sx(j)+tstep*real(dSx);
    Sy(j+1)=Sy(j)+tstep*real(dSy);
    Sz(j+1)=Sz(j)+tstep*real(dSz);

    Svector=[Sx(j+1),Sy(j+1),Sz(j+1)];
    Slength(j)=norm(Svector);
    SvectorNorm=Svector/norm(Svector);
    Sx(j+1)=SvectorNorm(1);
    Sy(j+1)=SvectorNorm(2);
    Sz(j+1)=SvectorNorm(3);
end