%The full Green's function which both contains the bare Green's function
%and the connection to the spin

g0less0(k,:)=(1i./(2*pi)).*fermi(k,:).*(g0(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+gS(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));
g0great0(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(g0(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+gS(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));

G0less0(k,:)=g0less0(k,:)-J.*Sz(j).*(1i./(2*pi)).*fermi(k,:).*(E01(k,:)+F01(k,:)+E10(k,:)+F10(k,:));
G0great0(k,:)=g0great0(k,:)+J.*Sz(j).*(1i./(2*pi)).*(1-fermi(k,:)).*(E01(k,:)+F01(k,:)+E10(k,:)+F10(k,:));

G1xless0(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sx(j).*(E00(k,:)+F00(k,:))-1i.*Sy(j).*(E10(k,:)+F10(k,:))+1i.*Sy(j).*(E01(k,:)+F01(k,:))+1i.*Sx(j).*(E11(k,:)+F11(k,:)));
G1xgreat0(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sx(j).*(E00(k,:)+F00(k,:))-1i.*Sy(j).*(E10(k,:)+F10(k,:))+1i.*Sy(j).*(E01(k,:)+F01(k,:))+1i.*Sx(j).*(E11(k,:)+F11(k,:)));

G1yless0(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sy(j).*(E00(k,:)+F00(k,:))+1i.*Sx(j).*(E10(k,:)+F10(k,:))-1i.*Sx(j).*(E01(k,:)+F01(k,:))+1i.*Sy(j).*(E11(k,:)+F11(k,:)));
G1ygreat0(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sy(j).*(E00(k,:)+F00(k,:))+1i.*Sx(j).*(E10(k,:)+F10(k,:))-1i.*Sx(j).*(E01(k,:)+F01(k,:))+1i.*Sy(j).*(E11(k,:)+F11(k,:)));

g1less0(k,:)=(1i./(2*pi)).*fermi(k,:).*(gS(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+g0(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));
g1great0(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(gS(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+g0(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));

G1zless0(k,:)=g1less0(k,:)-J.*Sz(j).*(1i./(2*pi)).*fermi(k,:).*(E00(k,:)+F00(k,:));
G1zgreat0(k,:)=g1great0(k,:)+J.*Sz(j).*(1i./(2*pi)).*(1-fermi(k,:)).*(E00(k,:)+F00(k,:));