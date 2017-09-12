%Integrating over energies for the different Green's functions

G0less(i)=trapz(w,G0less0(1,:)+G0less0(2,:));
G0great(i)=trapz(w,G0great0(1,:)+G0great0(2,:));

G1xless(i)=trapz(w,G1xless0(1,:)+G1xless0(2,:));
G1xgreat(i)=trapz(w,G1xgreat0(1,:)+G1xgreat0(2,:));

G1yless(i)=trapz(w,G1yless0(1,:)+G1yless0(2,:));
G1ygreat(i)=trapz(w,G1ygreat0(1,:)+G1ygreat0(2,:));

G1zless(i)=trapz(w,G1zless0(1,:)+G1zless0(2,:));
G1zgreat(i)=trapz(w,G1zgreat0(1,:)+G1zgreat0(2,:));

Kless(i)=trapz(w,K0less(1,:)+K0less(2,:));
Kgreat(i)=trapz(w,K0great(1,:)+K0great(2,:));

G0less2(i)=trapz(w,G0less02(1,:)+G0less02(2,:));
G0great2(i)=trapz(w,G0great02(1,:)+G0great02(2,:));

G1xless2(i)=trapz(w,G1xless02(1,:)+G1xless02(2,:));
G1xgreat2(i)=trapz(w,G1xgreat02(1,:)+G1xgreat02(2,:));

G1yless2(i)=trapz(w,G1yless02(1,:)+G1yless02(2,:));
G1ygreat2(i)=trapz(w,G1ygreat02(1,:)+G1ygreat02(2,:));

G1zless2(i)=trapz(w,G1zless02(1,:)+G1zless02(2,:));
G1zgreat2(i)=trapz(w,G1zgreat02(1,:)+G1zgreat02(2,:));