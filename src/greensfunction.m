function [G0less, G0great, G1xless, G1xgreat, G1yless, G1ygreat, G1zless, G1zgreat] = greensfunction(t, tau, t0, t1, eps, gamma, gamma0, gammaS, eV, w, fermi, S, J)

%Defining some variables for speed
A=zeros(2,length(w));
B=zeros(2,length(w));
C=zeros(2,2,length(w));
D=zeros(2,2,length(w));
A0=zeros(2,length(w));
A1=zeros(2,length(w));
B0=zeros(2,length(w));
B1=zeros(2,length(w));
C00=zeros(2,length(w));
C01=zeros(2,length(w));
C10=zeros(2,length(w));
C11=zeros(2,length(w));
D00=zeros(2,length(w));
D01=zeros(2,length(w));
D10=zeros(2,length(w));
D11=zeros(2,length(w));
E00=zeros(2,length(w));
E01=zeros(2,length(w));
E10=zeros(2,length(w));
E11=zeros(2,length(w));
F00=zeros(2,length(w));
F01=zeros(2,length(w));
F10=zeros(2,length(w));
F11=zeros(2,length(w));
g0less0=zeros(2,length(w));
g0great0=zeros(2,length(w));
G0less0=zeros(2,length(w));
G0great0=zeros(2,length(w));
g1less0=zeros(2,length(w));
g1great0=zeros(2,length(w));
G1xless0=zeros(2,length(w));
G1xgreat0=zeros(2,length(w));
G1yless0=zeros(2,length(w));
G1ygreat0=zeros(2,length(w));
G1zless0=zeros(2,length(w));
G1zgreat0=zeros(2,length(w));

%Green's function for a voltage pulse
for k=1:2
    for m=1:2
        if tau < t0
                A(m,:) = -(1./2).*(exp(-1i.*w.*tau)./(w-eps(m)+1i.*gamma(m)/2));
        elseif tau < t1
                A(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*gamma(m)/2).*(tau-t0)-1i.*w.*t0).*(1./(w-eps(m)+1i.*gamma(m)/2)-1./(w+eV(k)-eps(m)+1i*gamma(m)/2))...
                +exp(-1i.*eV(k).*(tau-t0)-1i.*w.*tau)./(w+eV(k)-eps(m)+1i*gamma(m)/2));
        else
                A(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*gamma(m)/2).*(tau)).*(1./(w-eps(m)+1i.*gamma(m)/2).*(...
                exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t0)+exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*tau)-exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t1))...
                +1./(w+eV(k)-eps(m)+1i*gamma(m)/2).*(exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t0))));
        end

        if t < t0
                B(m,:) = -(1./2).*(exp(1i.*w.*t)./(w-eps(m)-1i.*gamma(m)/2));
        elseif t < t1
                B(m,:) = -(1./2).*(exp(-1i.*(eps(m)+1i.*gamma(m)/2).*(t0-t)+1i.*w.*t0).*(1./(w-eps(m)-1i.*gamma(m)/2)-1./(w+eV(k)-eps(m)-1i*gamma(m)/2))...
                +exp(1i.*eV(k).*(t-t0)+1i.*w.*t)./(w+eV(k)-eps(m)-1i*gamma(m)/2));
        else
                B(m,:) = -(1./2).*(exp(1i.*(eps(m)+1i.*gamma(m)/2).*(t)).*(1./(w-eps(m)-1i.*gamma(m)/2).*(...
                exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t0)+exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t)-exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t1))...
                +1./(w+eV(k)-eps(m)-1i*gamma(m)/2).*(exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t0))));
        end

        for l=1:2

            if tau < t0
                    C(m,l,:)=(1./4).*exp(-1i.*w.*tau)./((w-eps(m)+1i*gamma(m)/2).*(w-eps(l)+1i*gamma(l)/2));
            elseif tau < t1
                    C(m,l,:)=(1./4).*(exp(-1i.*w.*t0-1i.*(eps(m)-1i.*gamma(m)/2).*(tau-t0))./((w-eps(m)+1i*gamma(m)/2).*(w-eps(l)+1i*gamma(l)/2))...
                    +(exp(-1i.*w.*tau-1i.*eV(k).*(tau-t0))-exp(-1i.*w.*t0-1i.*(eps(m)-1i.*gamma(m)/2).*(tau-t0)))./((w+eV(k)-eps(m)+1i*gamma(m)/2).*(w+eV(k)-eps(l)+1i*gamma(l)/2)));
            else
                    C(m,l,:) = (1./4).*(exp(-1i.*(eps(m)-1i.*gamma(m)/2).*(tau)).*(1./((w-eps(m)+1i*gamma(m)/2).*(w-eps(l)+1i*gamma(l)/2)).*(...
                    exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t0)+exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*tau)-exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t1))...
                    +1./((w+eV(k)-eps(m)+1i*gamma(m)/2).*(w+eV(k)-eps(l)+1i*gamma(l)/2)).*(exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*gamma(m)/2-w).*t0))));
            end

            if t < t0
                    D(m,l,:)=(1./4).*(exp(1i.*w.*t)./((w-eps(m)-1i*gamma(m)/2).*(w-eps(l)-1i*gamma(l)/2)));
            elseif t < t1
                    D(m,l,:)=(1./4).*(exp(1i.*w.*t0-1i.*(eps(m)+1i.*gamma(m)/2).*(t0-t))./((w-eps(m)-1i*gamma(m)/2).*(w-eps(l)-1i*gamma(l)/2))...
                    +(exp(1i.*w.*t+1i.*eV(k).*(t-t0))-exp(1i.*w.*t0-1i.*(eps(m)+1i.*gamma(m)/2).*(t0-t)))./((w+eV(k)-eps(m)-1i*gamma(m)/2).*(w+eV(k)-eps(l)-1i*gamma(l)/2)));
            else
                    D(m,l,:) = (1./4).*(exp(1i.*(eps(m)+1i.*gamma(m)/2).*(t)).*(1./((w-eps(m)-1i*gamma(m)/2).*(w-eps(l)-1i*gamma(l)/2)).*(...
                    exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t0)+exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t)-exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t1))...
                    +1./((w+eV(k)-eps(m)-1i*gamma(m)/2).*(w+eV(k)-eps(l)-1i*gamma(l)/2)).*(exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*gamma(m)/2-w).*t0))));
            end
        end
    end

    A0(k,:)=A(1,:)+A(2,:);
    A1(k,:)=A(1,:)-A(2,:);

    B0(k,:)=B(1,:)+B(2,:);
    B1(k,:)=B(1,:)-B(2,:);

    C00(k,:)=C(1,1,:)+C(1,2,:)+C(2,1,:)+C(2,2,:);
    C01(k,:)=C(1,1,:)-C(1,2,:)+C(2,1,:)-C(2,2,:);
    C10(k,:)=C(1,1,:)+C(1,2,:)-C(2,1,:)-C(2,2,:);
    C11(k,:)=C(1,1,:)-C(1,2,:)-C(2,1,:)+C(2,2,:);

    D00(k,:)=D(1,1,:)+D(1,2,:)+D(2,1,:)+D(2,2,:);
    D01(k,:)=D(1,1,:)-D(1,2,:)+D(2,1,:)-D(2,2,:);
    D10(k,:)=D(1,1,:)+D(1,2,:)-D(2,1,:)-D(2,2,:);
    D11(k,:)=D(1,1,:)-D(1,2,:)-D(2,1,:)+D(2,2,:);

    E00(k,:)=gamma0(k).*(C00(k,:).*B0(k,:)+C01(k,:).*B1(k,:))+gammaS(k).*(C00(k,:).*B1(k,:)+C01(k,:).*B0(k,:));
    E01(k,:)=gammaS(k).*(C00(k,:).*B0(k,:)+C01(k,:).*B1(k,:))+gamma0(k).*(C00(k,:).*B1(k,:)+C01(k,:).*B0(k,:));
    E10(k,:)=gamma0(k).*(C10(k,:).*B0(k,:)+C11(k,:).*B1(k,:))+gammaS(k).*(C10(k,:).*B1(k,:)+C11(k,:).*B0(k,:));
    E11(k,:)=gammaS(k).*(C10(k,:).*B0(k,:)+C11(k,:).*B1(k,:))+gamma0(k).*(C10(k,:).*B1(k,:)+C11(k,:).*B0(k,:));
    F00(k,:)=gamma0(k).*(A0(k,:).*D00(k,:)+A1(k,:).*D01(k,:))+gammaS(k).*(A0(k,:).*D01(k,:)+A1(k,:).*D00(k,:));
    F01(k,:)=gamma0(k).*(A0(k,:).*D10(k,:)+A1(k,:).*D11(k,:))+gammaS(k).*(A0(k,:).*D11(k,:)+A1(k,:).*D10(k,:));
    F10(k,:)=gammaS(k).*(A0(k,:).*D00(k,:)+A1(k,:).*D01(k,:))+gamma0(k).*(A0(k,:).*D01(k,:)+A1(k,:).*D00(k,:));
    F11(k,:)=gammaS(k).*(A0(k,:).*D10(k,:)+A1(k,:).*D11(k,:))+gamma0(k).*(A0(k,:).*D11(k,:)+A1(k,:).*D10(k,:));

    %The full Green's function which both contains the bare Green's function
    %and the connection to the spin

    g0less0(k,:)=(1i./(2*pi)).*fermi(k,:).*(gamma0(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+gammaS(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));
    g0great0(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(gamma0(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+gammaS(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));

    G0less0(k,:)=g0less0(k,:)-J.*S(3).*(1i./(2*pi)).*fermi(k,:).*(E01(k,:)+F01(k,:)+E10(k,:)+F10(k,:));
    G0great0(k,:)=g0great0(k,:)+J.*S(3).*(1i./(2*pi)).*(1-fermi(k,:)).*(E01(k,:)+F01(k,:)+E10(k,:)+F10(k,:));

    G1xless0(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(S(1).*(E00(k,:)+F00(k,:))-1i.*S(2).*(E10(k,:)+F10(k,:))+1i.*S(2).*(E01(k,:)+F01(k,:))+1i.*S(1).*(E11(k,:)+F11(k,:)));
    G1xgreat0(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(S(1).*(E00(k,:)+F00(k,:))-1i.*S(2).*(E10(k,:)+F10(k,:))+1i.*S(2).*(E01(k,:)+F01(k,:))+1i.*S(1).*(E11(k,:)+F11(k,:)));

    G1yless0(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(S(2).*(E00(k,:)+F00(k,:))+1i.*S(1).*(E10(k,:)+F10(k,:))-1i.*S(1).*(E01(k,:)+F01(k,:))+1i.*S(2).*(E11(k,:)+F11(k,:)));
    G1ygreat0(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(S(2).*(E00(k,:)+F00(k,:))+1i.*S(1).*(E10(k,:)+F10(k,:))-1i.*S(1).*(E01(k,:)+F01(k,:))+1i.*S(2).*(E11(k,:)+F11(k,:)));

    g1less0(k,:)=(1i./(2*pi)).*fermi(k,:).*(gammaS(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+gamma0(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));
    g1great0(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(gammaS(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+gamma0(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));

    G1zless0(k,:)=g1less0(k,:)-J.*S(3).*(1i./(2*pi)).*fermi(k,:).*(E00(k,:)+F00(k,:)+E11(k,:)+F11(k,:));
    G1zgreat0(k,:)=g1great0(k,:)+J.*S(3).*(1i./(2*pi)).*(1-fermi(k,:)).*(E00(k,:)+F00(k,:)+E11(k,:)+F11(k,:));
end

%Integrating over energies for the different Green's functions

G0less=trapz(w,G0less0(1,:)+G0less0(2,:));
G0great=trapz(w,G0great0(1,:)+G0great0(2,:));

G1xless=trapz(w,G1xless0(1,:)+G1xless0(2,:));
G1xgreat=trapz(w,G1xgreat0(1,:)+G1xgreat0(2,:));

G1yless=trapz(w,G1yless0(1,:)+G1yless0(2,:));
G1ygreat=trapz(w,G1ygreat0(1,:)+G1ygreat0(2,:));

G1zless=trapz(w,G1zless0(1,:)+G1zless0(2,:));
G1zgreat=trapz(w,G1zgreat0(1,:)+G1zgreat0(2,:));
end
