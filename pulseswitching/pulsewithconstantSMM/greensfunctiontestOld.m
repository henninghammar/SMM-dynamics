% Test of the greens function
% ================================================
% Author: Henning Hammar
% ------------------

clear

addpath('../../src')

tic

%Initial values
pL=0; %Polarization of gamma_up and gamma_down
pR=0;
gamma=1;
g0(1)=gamma;
g0(2)=gamma;
eV(1)=gamma; %bias voltage on left lead
eV(2)=-gamma; %bias voltage on right lead
gfactor=2; %g-factor
myB=5.78838175*10^(-2); %Bohr magneton in meV*T^-1
B=1; %Magnetic field in Tesla
J=0.1*gamma; %Coupling strength
Sx = 0;
Sy = sin(pi/4);
Sz = cos(pi/4);
wL=gfactor*myB*B; %Frequency
epsilon=0; %Energy level of the quantum dot
eps(1)=epsilon+0.5*wL;
eps(2)=epsilon-0.5*wL;

%Time variables and time and energy step-size
tstep=0.1;
tstep2=0.1/gamma;
tback=200*gamma;
t0=0;%50*tstep;
t1=3;

kB=8.617324*10^-2; %Boltzmanns constant, in meV*K^-1
T(1)=1; %Temperature in K
T(2)=1;

%Integration values in time and energy
t=5;
step=0.1;
w=[-20:step:20];

%Initial values
gS(1)=pL*g0(1);
gS(2)=pR*g0(2);
g0tot=g0(1)+g0(2); %Gamma0 in eV (Coupling strength between leads and the dot)
gStot=gS(1)+gS(2);
g(1)=g0(1)/2*(1+pL)+g0(2)/2*(1+pR);
g(2)=g0(1)/2*(1-pL)+g0(2)/2*(1-pR);

step=0.1;
w=[-20:step:20];

fermi=zeros(2,length(w));
beta(1)=1/(kB*T(1)); %Beta value
beta(2)=1/(kB*T(2)); %Beta value
fermi(1,:)=1./(1+exp(beta(1).*w));
fermi(2,:)=1./(1+exp(beta(2).*w));

tau = [-tstep2*tback+t:tstep2:t];

%Defining some variables
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
G0less=zeros(1,tback+1);
G0great=zeros(1,tback+1);
G1xless=zeros(1,tback+1);
G1xgreat=zeros(1,tback+1);
G1yless=zeros(1,tback+1);
G1ygreat=zeros(1,tback+1);
G1zless=zeros(1,tback+1);
G1zgreat=zeros(1,tback+1);
A2=zeros(2,length(w));
B2=zeros(2,length(w));
C2=zeros(2,2,length(w));
D2=zeros(2,2,length(w));
A02=zeros(2,length(w));
A12=zeros(2,length(w));
B02=zeros(2,length(w));
B12=zeros(2,length(w));
C002=zeros(2,length(w));
C012=zeros(2,length(w));
C102=zeros(2,length(w));
C112=zeros(2,length(w));
D002=zeros(2,length(w));
D012=zeros(2,length(w));
D102=zeros(2,length(w));
D112=zeros(2,length(w));
E002=zeros(2,length(w));
E012=zeros(2,length(w));
E102=zeros(2,length(w));
E112=zeros(2,length(w));
F002=zeros(2,length(w));
F012=zeros(2,length(w));
F102=zeros(2,length(w));
F112=zeros(2,length(w));
g0less02=zeros(2,length(w));
g0great02=zeros(2,length(w));
G0less02=zeros(2,length(w));
G0great02=zeros(2,length(w));
g1less02=zeros(2,length(w));
g1great02=zeros(2,length(w));
G1xless02=zeros(2,length(w));
G1xgreat02=zeros(2,length(w));
G1yless02=zeros(2,length(w));
G1ygreat02=zeros(2,length(w));
G1zless02=zeros(2,length(w));
G1zgreat02=zeros(2,length(w));
G0less2=zeros(1,tback+1);
G0great2=zeros(1,tback+1);
G1xless2=zeros(1,tback+1);
G1xgreat2=zeros(1,tback+1);
G1yless2=zeros(1,tback+1);
G1ygreat2=zeros(1,tback+1);
G1zless2=zeros(1,tback+1);
G1zgreat2=zeros(1,tback+1);
K0less=zeros(2,length(w));
K0great=zeros(2,length(w));
Kless=zeros(1,tback+1);
Kgreat=zeros(1,tback+1);
Ic0=zeros(1,length(t));
Ic1=zeros(1,length(t));
Isx0=zeros(1,length(t));
Isx1=zeros(1,length(t));
Isy0=zeros(1,length(t));
Isy1=zeros(1,length(t));
Isz0=zeros(1,length(t));
Isz1=zeros(1,length(t));

j = 1;

for i=1:tback+1

    for k=1:2

        %Function of (t',t)

        for m=1:2
            if tau(i) < t0
                    A(m,:) = -(1./2).*(exp(-1i.*w.*tau(i))./(w-eps(m)+1i.*g(m)));
            elseif tau(i) < t1
                    A(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*g(m)).*(tau(i)-t0)-1i.*w.*t0).*(1./(w-eps(m)+1i.*g(m))-1./(w+eV(k)-eps(m)+1i*g(m)))...
                    +exp(-1i.*eV(k).*(tau(i)-t0)-1i.*w.*tau(i))./(w+eV(k)-eps(m)+1i*g(m)));
            else
                    A(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*g(m)).*(tau(i))).*(1./(w-eps(m)+1i.*g(m)).*(...
                    exp(1i.*(eps(m)-1i.*g(m)-w).*t0)+exp(1i.*(eps(m)-1i.*g(m)-w).*tau(i))-exp(1i.*(eps(m)-1i.*g(m)-w).*t1))...
                    +1./(w+eV(k)-eps(m)+1i*g(m)).*(exp(1i.*(eps(m)-1i.*g(m)-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*g(m)-w).*t0))));
            end

            if t < t0
                    B(m,:) = -(1./2).*(exp(1i.*w.*t)./(w-eps(m)-1i.*g(m)));
            elseif t < t1
                    B(m,:) = -(1./2).*(exp(-1i.*(eps(m)+1i.*g(m)).*(t0-t)+1i.*w.*t0).*(1./(w-eps(m)-1i.*g(m))-1./(w+eV(k)-eps(m)-1i*g(m)))...
                    +exp(1i.*eV(k).*(t-t0)+1i.*w.*t)./(w+eV(k)-eps(m)-1i*g(m)));
            else
                    B(m,:) = -(1./2).*(exp(1i.*(eps(m)+1i.*g(m)).*(t)).*(1./(w-eps(m)-1i.*g(m)).*(...
                    exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*t)-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
                    +1./(w+eV(k)-eps(m)-1i*g(m)).*(exp(-1i.*(eps(m)+1i.*g(m)-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t0))));
            end

            for l=1:2

                if tau(i) < t0
                        C(m,l,:)=(1./4).*exp(-1i.*w.*tau(i))./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l)));
                elseif tau(i) < t1
                        C(m,l,:)=(1./4).*(exp(-1i.*w.*t0-1i.*(eps(m)-1i.*g(m)).*(tau(i)-t0))./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l)))...
                        +(exp(-1i.*w.*tau(i)-1i.*eV(k).*(tau(i)-t0))-exp(-1i.*w.*t0-1i.*(eps(m)-1i.*g(m)).*(tau(i)-t0)))./((w+eV(k)-eps(m)+1i*g(m)).*(w+eV(k)-eps(l)+1i*g(l))));
                else
                        C(m,l,:) = (1./4).*(exp(-1i.*(eps(m)-1i.*g(m)).*(tau(i))).*(1./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l))).*(...
                        exp(1i.*(eps(m)-1i.*g(m)-w).*t0)+exp(1i.*(eps(m)-1i.*g(m)-w).*tau(i))-exp(1i.*(eps(m)-1i.*g(m)-w).*t1))...
                        +1./((w+eV(k)-eps(m)+1i*g(m)).*(w+eV(k)-eps(l)+1i*g(l))).*(exp(1i.*(eps(m)-1i.*g(m)-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*g(m)-w).*t0))));
                end

                if t < t0
                        D(m,l,:)=(1./4).*(exp(1i.*w.*t)./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))));
                elseif t < t1
                        D(m,l,:)=(1./4).*(exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-t))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l)))...
                        +(exp(1i.*w.*t+1i.*eV(k).*(t-t0))-exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-t)))./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))));
                else
                        D(m,l,:) = (1./4).*(exp(1i.*(eps(m)+1i.*g(m)).*(t)).*(1./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))).*(...
                        exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*t)-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
                        +1./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))).*(exp(-1i.*(eps(m)+1i.*g(m)-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t0))));
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

        E00(k,:)=g0(k).*(C00(k,:).*B0(k,:)+C01(k,:).*B1(k,:))+gS(k).*(C00(k,:).*B1(k,:)+C01(k,:).*B0(k,:));
        E01(k,:)=gS(k).*(C00(k,:).*B0(k,:)+C01(k,:).*B1(k,:))+g0(k).*(C00(k,:).*B1(k,:)+C01(k,:).*B0(k,:));
        E10(k,:)=g0(k).*(C10(k,:).*B0(k,:)+C11(k,:).*B1(k,:))+gS(k).*(C10(k,:).*B1(k,:)+C11(k,:).*B0(k,:));
        E11(k,:)=gS(k).*(C10(k,:).*B0(k,:)+C11(k,:).*B1(k,:))+g0(k).*(C10(k,:).*B1(k,:)+C11(k,:).*B0(k,:));
        F00(k,:)=g0(k).*(A0(k,:).*D00(k,:)+A1(k,:).*D01(k,:))+gS(k).*(A0(k,:).*D01(k,:)+A1(k,:).*D00(k,:));
        F01(k,:)=g0(k).*(A0(k,:).*D10(k,:)+A1(k,:).*D11(k,:))+gS(k).*(A0(k,:).*D11(k,:)+A1(k,:).*D10(k,:));
        F10(k,:)=gS(k).*(A0(k,:).*D00(k,:)+A1(k,:).*D01(k,:))+g0(k).*(A0(k,:).*D01(k,:)+A1(k,:).*D00(k,:));
        F11(k,:)=gS(k).*(A0(k,:).*D10(k,:)+A1(k,:).*D11(k,:))+g0(k).*(A0(k,:).*D11(k,:)+A1(k,:).*D10(k,:));


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
    end

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
end

%Savedata
outputFolder = 'output';
outputFilename = sprintf('%s/test10.mat', outputFolder);
%save(outputFilename)
