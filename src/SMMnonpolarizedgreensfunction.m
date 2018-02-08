function [G0less, G0great, G1xless, G1xgreat, G1yless, G1ygreat, G1zless, G1zgreat] = SMMnonpolarizedgreensfunction(t, tau, t0, t1, epsilon0, gamma0, gamma0lead, mu, w, fermi, S, J)

%Defining some variables for speed
A=zeros(2,length(w));
B=zeros(2,length(w));
C=zeros(2,length(w));
D=zeros(2,length(w));
G0less0=zeros(2,length(w));
G0great0=zeros(2,length(w));
G1lesskernel=zeros(2,length(w));
G1greatkernel=zeros(2,length(w));
G1xless0=zeros(2,length(w));
G1xgreat0=zeros(2,length(w));
G1yless0=zeros(2,length(w));
G1ygreat0=zeros(2,length(w));
G1zless0=zeros(2,length(w));
G1zgreat0=zeros(2,length(w));

%Green's function for a voltage pulse
for k=1:2
    if tau < t0
            A(k,:) = -exp(-1i.*w.*tau)./(w-epsilon0+1i.*gamma0/4);
    elseif tau < t1
            A(k,:) = -(exp(-1i.*(epsilon0-1i.*gamma0/4).*(tau-t0)-1i.*w.*t0).*(...
            1./(w-epsilon0+1i.*gamma0/4)-1./(w+mu(k)-epsilon0+1i*gamma0/4))...
            +exp(-1i.*mu(k).*(tau-t0)-1i.*w.*tau)./(w+mu(k)-epsilon0+1i*gamma0/4));
    else
            A(k,:) = -exp(-1i.*(epsilon0-1i.*gamma0/4).*tau).*(...
            1./(w-epsilon0+1i.*gamma0/4).*(...
            exp(1i.*(epsilon0-1i.*gamma0/4-w).*t0)...
            +exp(1i.*(epsilon0-1i.*gamma0/4-w).*tau)...
            -exp(1i.*(epsilon0-1i.*gamma0/4-w).*t1))...
            +1./(w+mu(k)-epsilon0+1i*gamma0/4).*(...
            exp(1i.*(epsilon0-1i.*gamma0/4-w).*t1-1i.*mu(k).*(t1-t0))...
            -exp(1i.*(epsilon0-1i.*gamma0/4-w).*t0)));
    end

    if t < t0
            B(k,:) = -exp(1i.*w.*t)./(w-epsilon0-1i.*gamma0/4);
    elseif t < t1
            B(k,:) = -(exp(-1i.*(epsilon0+1i.*gamma0/4).*(t0-t)+1i.*w.*t0).*(...
            1./(w-epsilon0-1i.*gamma0/4)-1./(w+mu(k)-epsilon0-1i*gamma0/4))...
            +exp(1i.*mu(k).*(t-t0)+1i.*w.*t)./(w+mu(k)-epsilon0-1i*gamma0/4));
    else
            B(k,:) = -exp(1i.*(epsilon0+1i.*gamma0/4).*t).*(...
            1./(w-epsilon0-1i.*gamma0/4).*(...
            exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t0)...
            +exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t)...
            -exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t1))...
            +1./(w+mu(k)-epsilon0-1i*gamma0/4).*(...
            exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t1+1i.*mu(k).*(t1-t0))...
            -exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t0)));
    end

    if tau < t0
            C(k,:)=exp(-1i.*w.*tau)./((w-epsilon0+1i*gamma0/4).*(w-epsilon0+1i*gamma0/4));
    elseif tau < t1
            C(k,:)=exp(-1i.*w.*t0-1i.*(epsilon0-1i.*gamma0/4).*(tau-t0))./((w-epsilon0+1i*gamma0/4).*(w-epsilon0+1i*gamma0/4))...
            +(exp(-1i.*w.*tau-1i.*mu(k).*(tau-t0))...
            -exp(-1i.*w.*t0-1i.*(epsilon0-1i.*gamma0/4).*(tau-t0)))./((w+mu(k)-epsilon0+1i*gamma0/4).*(w+mu(k)-epsilon0+1i*gamma0/4));
    else
            C(k,:) = exp(-1i.*(epsilon0-1i.*gamma0/4).*tau).*(...
            1./((w-epsilon0+1i*gamma0/4).*(w-epsilon0+1i*gamma0/4)).*(...
            exp(1i.*(epsilon0-1i.*gamma0/4-w).*t0)...
            +exp(1i.*(epsilon0-1i.*gamma0/4-w).*tau)...
            -exp(1i.*(epsilon0-1i.*gamma0/4-w).*t1))...
            +1./((w+mu(k)-epsilon0+1i*gamma0/4).*(w+mu(k)-epsilon0+1i*gamma0/4)).*(...
            exp(1i.*(epsilon0-1i.*gamma0/4-w).*t1-1i.*mu(k).*(t1-t0))...
            -exp(1i.*(epsilon0-1i.*gamma0/4-w).*t0)));
    end

    if t < t0
            D(k,:)=exp(1i.*w.*t)./((w-epsilon0-1i*gamma0/4).*(w-epsilon0-1i*gamma0/4));
    elseif t < t1
            D(k,:)=exp(1i.*w.*t0+1i.*(epsilon0+1i.*gamma0/4).*(t-t0))./((w-epsilon0-1i*gamma0/4).*(w-epsilon0-1i*gamma0/4))...
            +(exp(1i.*w.*t+1i.*mu(k).*(t-t0))...
            -exp(1i.*w.*t0+1i.*(epsilon0+1i.*gamma0/4).*(t-t0)))./((w+mu(k)-epsilon0-1i*gamma0/4).*(w+mu(k)-epsilon0-1i*gamma0/4));
    else
            D(k,:) = exp(1i.*(epsilon0+1i.*gamma0/4).*(t)).*(...
            1./((w-epsilon0-1i*gamma0/4).*(w-epsilon0-1i*gamma0/4)).*(...
            exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t0)...
            +exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t)...
            -exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t1))...
            +1./((w+mu(k)-epsilon0-1i*gamma0/4).*(w+mu(k)-epsilon0-1i*gamma0/4)).*(...
            exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t1+1i.*mu(k).*(t1-t0))...
            -exp(-1i.*(epsilon0+1i.*gamma0/4-w).*t0)));
    end

    %The full Green's function which both contains the bare Green's function
    %and the connection to the spin
    G0less0(k,:)=(1i./(2*pi)).*fermi(k,:).*(gamma0lead(k).*A(k,:).*B(k,:));
    G0great0(k,:)=(-1i./(2*pi)).*(1-fermi(k,:)).*(gamma0lead(k).*A(k,:).*B(k,:));

    G1lesskernel(k,:) = -(1i./(2*pi)).*J.*fermi(k,:).*gamma0lead(k).*(C(k,:).*B(k,:)+A(k,:).*D(k,:));
    G1greatkernel(k,:) = (1i./(2*pi)).*J.*(1-fermi(k,:)).*gamma0lead(k).*(C(k,:).*B(k,:)+A(k,:).*D(k,:));

    G1xless0(k,:) = S(1)*G1lesskernel(k,:);
    G1xgreat0(k,:) = S(1)*G1greatkernel(k,:);

    G1yless0(k,:) = S(2)*G1lesskernel(k,:);
    G1ygreat0(k,:) = S(2)*G1greatkernel(k,:);

    G1zless0(k,:) = S(3)*G1lesskernel(k,:);
    G1zgreat0(k,:) = S(3)*G1greatkernel(k,:);
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
