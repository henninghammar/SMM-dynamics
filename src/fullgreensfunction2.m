%The full Green's function which both contains the bare Green's function
%and the connection to the spin

g0less02(k,:)=(1i./(2*pi)).*fermi(k,:).*(g0(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+gS(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));
g0great02(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(g0(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+gS(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));

G0less02(k,:)=g0less02(k,:)-J.*Sz(j).*(1i./(2*pi)).*fermi(k,:).*(E012(k,:)+F012(k,:)+E102(k,:)+F102(k,:));
G0great02(k,:)=g0great02(k,:)+J.*Sz(j).*(1i./(2*pi)).*(1-fermi(k,:)).*(E012(k,:)+F012(k,:)+E102(k,:)+F102(k,:));

G1xless02(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sx(j).*(E002(k,:)+F002(k,:))-1i.*Sy(j).*(E102(k,:)+F102(k,:))+1i.*Sy(j).*(E012(k,:)+F012(k,:))+1i.*Sx(j).*(E112(k,:)+F112(k,:)));
G1xgreat02(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sx(j).*(E002(k,:)+F002(k,:))-1i.*Sy(j).*(E102(k,:)+F102(k,:))+1i.*Sy(j).*(E012(k,:)+F012(k,:))+1i.*Sx(j).*(E112(k,:)+F112(k,:)));

G1yless02(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sy(j).*(E002(k,:)+F002(k,:))+1i.*Sx(j).*(E102(k,:)+F102(k,:))-1i.*Sx(j).*(E012(k,:)+F012(k,:))+1i.*Sy(j).*(E112(k,:)+F112(k,:)));
G1ygreat02(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sy(j).*(E002(k,:)+F002(k,:))+1i.*Sx(j).*(E102(k,:)+F102(k,:))-1i.*Sx(j).*(E012(k,:)+F012(k,:))+1i.*Sy(j).*(E112(k,:)+F112(k,:)));

g1less02(k,:)=(1i./(2*pi)).*fermi(k,:).*(gS(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+g0(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));
g1great02(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(gS(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+g0(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));

G1zless02(k,:)=g1less02(k,:)-J.*Sz(j).*(1i./(2*pi)).*fermi(k,:).*(E002(k,:)+F002(k,:));
G1zgreat02(k,:)=g1great02(k,:)+J.*Sz(j).*(1i./(2*pi)).*(1-fermi(k,:)).*(E002(k,:)+F002(k,:));