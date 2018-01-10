%The core program

%Initial values
gS(1)=pL*g0(1);
gS(2)=pR*g0(2);
g0tot=g0(1)+g0(2); %Gamma0 in eV (Coupling strength between leads and the dot)
gStot=gS(1)+gS(2);
g(1)=g0(1)/2*(1+pL)+g0(2)/2*(1+pR);
g(2)=g0(1)/2*(1-pL)+g0(2)/2*(1-pR);

tstep2=0.1;
%tneg=[-tstep*tback-t0:tstep2:tmax];
t=[0:tstep:tmax];
step=0.1;
w=[-20:step:20];

fermi=zeros(2,length(w));
beta(1)=1/(kB*T(1)); %Beta value
beta(2)=1/(kB*T(2)); %Beta value
fermi(1,:)=1./(1+exp(beta(1).*w));
fermi(2,:)=1./(1+exp(beta(2).*w));

%Spin
Sx=zeros(1,length(t));
Sy=zeros(1,length(t));
Sz=zeros(1,length(t));
Sx(1)=-Sxy*sin(wL*t0);
Sy(1)=Sxy*cos(wL*t0);
Sz(1)=Sz0;

SxNonorm=zeros(1,length(t));
SyNonorm=zeros(1,length(t));
SzNonorm=zeros(1,length(t));
SxNonorm(1)=-Sxy*sin(wL*t0);
SyNonorm(1)=Sxy*cos(wL*t0);
SzNonorm(1)=Sz0;
dSxvect1=zeros(1,length(t));
dSyvect1=zeros(1,length(t));
dSzvect1=zeros(1,length(t));
dSxvect2=zeros(1,length(t));
dSyvect2=zeros(1,length(t));
dSzvect2=zeros(1,length(t));
dSxvectTot=zeros(1,length(t));
dSyvectTot=zeros(1,length(t));
dSzvectTot=zeros(1,length(t));

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

for j=1:length(t)

    if t(j)>=0
        J=J0;
    else
        J=0;
    end

    if(t(j)==floor(t(j)))
        timestep=t(j);
        timeused=toc;
        disp(['Timestep: ' num2str(timestep)])
        disp(['Timeused: ' num2str(timeused)])
    end

    %tau=tneg(j:j+tback);
    tau = [-tstep2*tback+t(j):tstep2:t(j)];

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

                if t(j) < t0
                        B(m,:) = -(1./2).*(exp(1i.*w.*t(j))./(w-eps(m)-1i.*g(m)));
                elseif t(j) < t1
                        B(m,:) = -(1./2).*(exp(-1i.*(eps(m)+1i.*g(m)).*(t0-t(j))+1i.*w.*t0).*(1./(w-eps(m)-1i.*g(m))-1./(w+eV(k)-eps(m)-1i*g(m)))...
                        +exp(1i.*eV(k).*(t(j)-t0)+1i.*w.*t(j))./(w+eV(k)-eps(m)-1i*g(m)));
                else
                        B(m,:) = -(1./2).*(exp(1i.*(eps(m)+1i.*g(m)).*(t(j))).*(1./(w-eps(m)-1i.*g(m)).*(...
                        exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*t(j))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
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

                    if t(j) < t0
                            D(m,l,:)=(1./4).*(exp(1i.*w.*t(j))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))));
                    elseif t(j) < t1
                            D(m,l,:)=(1./4).*(exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-t(j)))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l)))...
                            +(exp(1i.*w.*t(j)+1i.*eV(k).*(t(j)-t0))-exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-t(j))))./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))));
                    else
                            D(m,l,:) = (1./4).*(exp(1i.*(eps(m)+1i.*g(m)).*(t(j))).*(1./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))).*(...
                            exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*t(j))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
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

                if t(j) < t0
                    if tau(i) < t0
                        K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                        K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                    else
                        K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t0-tau(i)));
                        K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t0-tau(i)));
                    end
                elseif t(j) < t1
                    if tau(i) < t0
                        K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t0));
                        K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t0));
                    elseif tau(i) < t1
                        K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-tau(i)));
                        K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-tau(i)));
                    else
                        K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t1));
                        K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t1));
                    end
                else
                    if tau(i) < t1
                        K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t1-tau(i)));
                        K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t1-tau(i)));
                    else
                        K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                        K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                    end
                end


%             %Function of (t,t')

            for m=1:2
                if t(j) < t0
                        A2(m,:) = -(1./2).*(exp(-1i.*w.*t(j))./(w-eps(m)+1i.*g(m)));
                elseif t(j) < t1
                        A2(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*g(m)).*(t(j)-t0)-1i.*w.*t0).*(1./(w-eps(m)+1i.*g(m))-1./(w+eV(k)-eps(m)+1i*g(m)))...
                        +exp(-1i.*eV(k).*(t(j)-t0)-1i.*w.*t(j))./(w+eV(k)-eps(m)+1i*g(m)));
                else
                        A2(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*g(m)).*(t(j))).*(1./(w-eps(m)+1i.*g(m)).*(...
                        exp(1i.*(eps(m)-1i.*g(m)-w).*t0)+exp(1i.*(eps(m)-1i.*g(m)-w).*t(j))-exp(1i.*(eps(m)-1i.*g(m)-w).*t1))...
                        +1./(w+eV(k)-eps(m)+1i*g(m)).*(exp(1i.*(eps(m)-1i.*g(m)-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*g(m)-w).*t0))));
                end

                if tau(i) < t0
                        B2(m,:) = -(1./2).*(exp(1i.*w.*tau(i))./(w-eps(m)-1i.*g(m)));
                elseif tau(i) < t1
                        B2(m,:) = -(1./2).*(exp(-1i.*(eps(m)+1i.*g(m)).*(t0-tau(i))+1i.*w.*t0).*(1./(w-eps(m)-1i.*g(m))-1./(w+eV(k)-eps(m)-1i*g(m)))...
                        +exp(1i.*eV(k).*(tau(i)-t0)+1i.*w.*tau(i))./(w+eV(k)-eps(m)-1i*g(m)));
                else
                        B2(m,:) = -(1./2).*(exp(1i.*(eps(m)+1i.*g(m)).*(tau(i))).*(1./(w-eps(m)-1i.*g(m)).*(...
                        exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*tau(i))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
                        +1./(w+eV(k)-eps(m)-1i*g(m)).*(exp(-1i.*(eps(m)+1i.*g(m)-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t0))));
                end

                for l=1:2

                    if t(j) < t0
                            C2(m,l,:)=(1./4).*exp(-1i.*w.*t(j))./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l)));
                    elseif t(j) < t1
                            C2(m,l,:)=(1./4).*(exp(-1i.*w.*t0-1i.*(eps(m)-1i.*g(m)).*(t(j)-t0))./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l)))...
                            +(exp(-1i.*w.*t(j)-1i.*eV(k).*(t(j)-t0))-exp(-1i.*w.*t0-1i.*(eps(m)-1i.*g(m)).*(t(j)-t0)))./((w+eV(k)-eps(m)+1i*g(m)).*(w+eV(k)-eps(l)+1i*g(l))));
                    else
                            C2(m,l,:) = (1./4).*(exp(-1i.*(eps(m)-1i.*g(m)).*(t(j))).*(1./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l))).*(...
                            exp(1i.*(eps(m)-1i.*g(m)-w).*t0)+exp(1i.*(eps(m)-1i.*g(m)-w).*t(j))-exp(1i.*(eps(m)-1i.*g(m)-w).*t1))...
                            +1./((w+eV(k)-eps(m)+1i*g(m)).*(w+eV(k)-eps(l)+1i*g(l))).*(exp(1i.*(eps(m)-1i.*g(m)-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*g(m)-w).*t0))));
                    end

                    if tau(i) < t0
                            D2(m,l,:)=(1./4).*(exp(1i.*w.*tau(i))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))));
                    elseif tau(i) < t1
                            D2(m,l,:)=(1./4).*(exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-tau(i)))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l)))...
                            +(exp(1i.*w.*tau(i)+1i.*eV(k).*(tau(i)-t0))-exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-tau(i))))./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))));
                    else
                            D2(m,l,:) = (1./4).*(exp(1i.*(eps(m)+1i.*g(m)).*(tau(i))).*(1./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))).*(...
                            exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*tau(i))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
                            +1./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))).*(exp(-1i.*(eps(m)+1i.*g(m)-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t0))));
                    end
                end
            end

            A02(k,:)=A2(1,:)+A2(2,:);
            A12(k,:)=A2(1,:)-A2(2,:);

            B02(k,:)=B2(1,:)+B2(2,:);
            B12(k,:)=B2(1,:)-B2(2,:);

            C002(k,:)=C2(1,1,:)+C2(1,2,:)+C2(2,1,:)+C2(2,2,:);
            C012(k,:)=C2(1,1,:)-C2(1,2,:)+C2(2,1,:)-C2(2,2,:);
            C102(k,:)=C2(1,1,:)+C2(1,2,:)-C2(2,1,:)-C2(2,2,:);
            C112(k,:)=C2(1,1,:)-C2(1,2,:)-C2(2,1,:)+C2(2,2,:);

            D002(k,:)=D2(1,1,:)+D2(1,2,:)+D2(2,1,:)+D2(2,2,:);
            D012(k,:)=D2(1,1,:)-D2(1,2,:)+D2(2,1,:)-D2(2,2,:);
            D102(k,:)=D2(1,1,:)+D2(1,2,:)-D2(2,1,:)-D2(2,2,:);
            D112(k,:)=D2(1,1,:)-D2(1,2,:)-D2(2,1,:)+D2(2,2,:);

            E002(k,:)=g0(k).*(C002(k,:).*B02(k,:)+C012(k,:).*B12(k,:))+gS(k).*(C002(k,:).*B12(k,:)+C012(k,:).*B02(k,:));
            E012(k,:)=gS(k).*(C002(k,:).*B02(k,:)+C012(k,:).*B12(k,:))+g0(k).*(C002(k,:).*B12(k,:)+C012(k,:).*B02(k,:));
            E102(k,:)=g0(k).*(C102(k,:).*B02(k,:)+C112(k,:).*B12(k,:))+gS(k).*(C102(k,:).*B12(k,:)+C112(k,:).*B02(k,:));
            E112(k,:)=gS(k).*(C102(k,:).*B02(k,:)+C112(k,:).*B12(k,:))+g0(k).*(C102(k,:).*B12(k,:)+C112(k,:).*B02(k,:));
            F002(k,:)=g0(k).*(A02(k,:).*D002(k,:)+A12(k,:).*D012(k,:))+gS(k).*(A02(k,:).*D012(k,:)+A12(k,:).*D002(k,:));
            F012(k,:)=g0(k).*(A02(k,:).*D102(k,:)+A12(k,:).*D112(k,:))+gS(k).*(A02(k,:).*D112(k,:)+A12(k,:).*D102(k,:));
            F102(k,:)=gS(k).*(A02(k,:).*D002(k,:)+A12(k,:).*D012(k,:))+g0(k).*(A02(k,:).*D012(k,:)+A12(k,:).*D002(k,:));
            F112(k,:)=gS(k).*(A02(k,:).*D102(k,:)+A12(k,:).*D112(k,:))+g0(k).*(A02(k,:).*D112(k,:)+A12(k,:).*D102(k,:));

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

            %No change in voltage and no spin dependency
%             greensfunctionAlt1
            %No change in voltage
            %greensfunctionAlt2
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

    G0greattot(j,:)=G0great;
    G1xgreattot(j,:)=G1xgreat;
    G1ygreattot(j,:)=G1ygreat;
    G1zgreattot(j,:)=G1zgreat;
%
%     G0greatint(j)=trapz(tau,G0great);
%     G1xgreatint(j)=trapz(tau,G1xgreat);
%     G1ygreatint(j)=trapz(tau,G1ygreat);
%     G1zgreatint(j)=trapz(tau,G1zgreat);

    G0greattot2(j,:)=G0great2;
    G1xgreattot2(j,:)=G1xgreat2;
    G1ygreattot2(j,:)=G1ygreat2;
    G1zgreattot2(j,:)=G1zgreat2;

     Ic0bare=Kgreat.*G0less+Kless.*G0great;
     Ic1bare=Kgreat.*G1zless+Kless.*G1zgreat;
     Ic0(j)=-4.*g0(1)*imag(trapz(tau,Ic0bare));
     Ic1(j)=-4.*gS(1)*imag(trapz(tau,Ic1bare));

     Isx0bare=Kgreat.*G1xless+Kless.*G1xgreat;
     Isx1bare=-1i.*(Kgreat.*G1yless+Kless.*G1ygreat);
     Isx0(j)=-4.*g0(1)*imag(trapz(tau,Isx0bare));
     Isx1(j)=-4.*gS(1)*imag(trapz(tau,Isx1bare));

     Isy0bare=Kgreat.*G1yless+Kless.*G1ygreat;
     Isy1bare=1i.*(Kgreat.*G1xless+Kless.*G1xgreat);
     Isy0(j)=-4.*g0(1)*imag(trapz(tau,Isy0bare));
     Isy1(j)=-4.*gS(1)*imag(trapz(tau,Isy1bare));

     Isz0bare=Kgreat.*G1zless+Kless.*G1zgreat;
     Isz1bare=Kgreat.*G0less+Kless.*G0great;
     Isz0(j)=-4.*g0(1)*imag(trapz(tau,Isz0bare));
     Isz1(j)=-4.*gS(1)*imag(trapz(tau,Isz1bare));

     Nbare = heaviside(t(j)-tau).*(G0less - G0great);
     N(j) = -2/pi.*imag(trapz(tau, Nbare));

     Mxbare = heaviside(t(j)-tau).*(G1xless - G1xgreat);
     Mx(j) = -1/pi.*imag(trapz(tau, Mxbare));

     Mybare = heaviside(t(j)-tau).*(G1yless - G1ygreat);
     My(j) = -1/pi.*imag(trapz(tau, Mybare));

     Mzbare = heaviside(t(j)-tau).*(G1zless - G1zgreat);
     Mz(j) = -1/pi.*imag(trapz(tau, Mzbare));

     if j == 1
        SxT=-Sxy*sin(wL*tau);
        SyT=Sxy*cos(wL*tau);
        SzT=Sz0*ones(1,length(tau));
     else
        SxT=[SxT(2:end),Sx(j)];
        SyT=[SyT(2:end),Sy(j)];
        SzT=[SzT(2:end),Sz(j)];
     end

    jH=1i.*J.^2.*(G0less.*G0great2-G0great.*G0less2-G1xless.*G1xgreat2+G1xgreat.*G1xless2...
        -G1yless.*G1ygreat2+G1ygreat.*G1yless2-G1zless.*G1zgreat2+G1zgreat.*G1zless2);
    jDMx=-J.^2.*(G0less.*G1xgreat2-G0great.*G1xless2-G1xless.*G0great2+G1xgreat.*G0less2);
    jDMy=-J.^2.*(G0less.*G1ygreat2-G0great.*G1yless2-G1yless.*G0great2+G1ygreat.*G0less2);
    jDMz=-J.^2.*(G0less.*G1zgreat2-G0great.*G1zless2-G1zless.*G0great2+G1zgreat.*G0less2);
    jIxx=2*1i.*J.^2.*(G1xless.*G1xgreat2-G1xgreat.*G1xless2);
    jIyy=2*1i.*J.^2.*(G1yless.*G1ygreat2-G1ygreat.*G1yless2);
    jIzz=2*1i.*J.^2.*(G1zless.*G1zgreat2-G1zgreat.*G1zless2);
    jIxy=1i.*J.^2.*(G1xless.*G1ygreat2-G1xgreat.*G1yless2+G1yless.*G1xgreat2-G1ygreat.*G1xless2);
    jIyx=jIxy;
    jIxz=1i.*J.^2.*(G1xless.*G1zgreat2-G1xgreat.*G1zless2+G1zless.*G1xgreat2-G1zgreat.*G1xless2);
    jIzx=jIxz;
    jIyz=1i.*J.^2.*(G1yless.*G1zgreat2-G1ygreat.*G1zless2+G1zless.*G1ygreat2-G1zgreat.*G1yless2);
    jIzy=jIyz;

    ejxbare=2*1i*J*((epsilon*G0less+wL/2.*G1zless).*G1xgreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1xless2...
        +(epsilon*G1xless-1i.*G1yless).*G0great2-(epsilon*G1xgreat-1i.*G1ygreat).*G0less2...
        -1i.*(epsilon*G1yless+1i.*G1xless).*G1zgreat2+1i.*(epsilon*G1zless+wL/2).*G1ygreat2...
        +1i.*(epsilon*G1ygreat+1i.*G1xgreat).*G1zless2-1i.*(epsilon*G1zgreat+wL/2).*G1yless2);
    ejybare=2*1i*J*((epsilon*G0less+wL/2.*G1zless).*G1ygreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1yless2...
        +(epsilon*G1yless+1i.*G1xless).*G0great2-(epsilon*G1ygreat+1i.*G1xgreat).*G0less2...
        -1i.*(epsilon*G1zless+wL/2).*G1xgreat2+1i.*(epsilon*G1xless-1i.*G1yless).*G1zgreat2...
        +1i.*(epsilon*G1zgreat+wL/2).*G1xless2-1i.*(epsilon*G1xgreat-1i.*G1ygreat).*G1zless2);
    ejzbare=2*1i*J*((epsilon*G0less+wL/2.*G1zless).*G1zgreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1zless2...
        +(epsilon*G1zless+wL/2).*G0great2-(epsilon*G1zgreat2+wL/2).*G0less2...
        -1i.*(epsilon*G1xless-1i.*G1yless).*G1ygreat2+1i.*(epsilon*G1yless+1i.*G1xless).*G1xgreat2...
        +1i.*(epsilon*G1xgreat-1i.*G1ygreat).*G1yless2-1i.*(epsilon*G1ygreat+1i.*G1xgreat).*G1xless2);

    ejx=trapz(tau,ejxbare);
    ejy=trapz(tau,ejybare);
    ejz=trapz(tau,ejzbare);

    dSxbare=jH.*(Sy(j).*SzT-Sz(j).*SyT)+Sy(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy(j).*(jDMx.*SyT-jDMy.*SxT)+Sz(j).*(jDMz.*SxT-jDMx.*SzT);
    dSybare=jH.*(Sz(j).*SxT-Sx(j).*SzT)+Sz(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz(j).*(jDMy.*SzT-jDMz.*SyT)+Sx(j).*(jDMx.*SyT-jDMy.*SxT);
    dSzbare=jH.*(Sx(j).*SyT-Sy(j).*SxT)+Sx(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx(j).*(jDMz.*SxT-jDMx.*SzT)+Sy(j).*(jDMy.*SzT-jDMz.*SyT);
    dSx1=-(wL+ejz).*Sy(j)+ejy.*Sz(j)+trapz(tau,dSxbare);
    dSy1=-ejx.*Sz(j)+(wL+ejz).*Sx(j)+trapz(tau,dSybare);
    dSz1=-ejy.*Sx(j)+ejx.*Sy(j)+trapz(tau,dSzbare);

    Sejx(j)=-ejz.*Sy(j)+ejy.*Sz(j);
    Sejy(j)=-ejx.*Sz(j)+ejz.*Sx(j);
    Sejz(j)=-ejy.*Sx(j)+ejx.*Sy(j);

    SBx(j)=-(wL+ejz).*Sy(j)+ejy.*Sz(j);
    SBy(j)=-ejx.*Sz(j)+(wL+ejz).*Sx(j);
    SBz(j)=-ejy.*Sx(j)+ejx.*Sy(j);

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

    Beffx=ejx;
    Beffy=ejy;
    Beffz=(wL+ejz);

    %The fields acting on the spin
    ejvect(:,j)=[ejx, ejy, ejz];
    
    Beffvect(:,j)=[Beffx, Beffy, Beffz];

    jHt(j)=trapz(tau,jH);
    jDMvect(:,j)=[trapz(tau,jDMx), trapz(tau,jDMy), trapz(tau,jDMz)];
    jIvect(:,j)=[trapz(tau,jIxx), trapz(tau,jIyy), trapz(tau,jIzz), trapz(tau,jIxy), trapz(tau,jIxz), trapz(tau,jIyz)];

    Sx2=Sx(j)+dSx1;
    Sy2=Sy(j)+dSy1;
    Sz2=Sz(j)+dSz1;
    SxT2=[SxT(2:end),Sx2];
    SyT2=[SyT(2:end),Sy2];
    SzT2=[SzT(2:end),Sz2];

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

                if t(j) < t0
                        B(m,:) = -(1./2).*(exp(1i.*w.*t(j))./(w-eps(m)-1i.*g(m)));
                elseif t(j) < t1
                        B(m,:) = -(1./2).*(exp(-1i.*(eps(m)+1i.*g(m)).*(t0-t(j))+1i.*w.*t0).*(1./(w-eps(m)-1i.*g(m))-1./(w+eV(k)-eps(m)-1i*g(m)))...
                        +exp(1i.*eV(k).*(t(j)-t0)+1i.*w.*t(j))./(w+eV(k)-eps(m)-1i*g(m)));
                else
                        B(m,:) = -(1./2).*(exp(1i.*(eps(m)+1i.*g(m)).*(t(j))).*(1./(w-eps(m)-1i.*g(m)).*(...
                        exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*t(j))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
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

                    if t(j) < t0
                            D(m,l,:)=(1./4).*(exp(1i.*w.*t(j))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))));
                    elseif t(j) < t1
                            D(m,l,:)=(1./4).*(exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-t(j)))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l)))...
                            +(exp(1i.*w.*t(j)+1i.*eV(k).*(t(j)-t0))-exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-t(j))))./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))));
                    else
                            D(m,l,:) = (1./4).*(exp(1i.*(eps(m)+1i.*g(m)).*(t(j))).*(1./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))).*(...
                            exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*t(j))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
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

            G0less0(k,:)=g0less0(k,:)-J.*Sz2.*(1i./(2*pi)).*fermi(k,:).*(E01(k,:)+F01(k,:)+E10(k,:)+F10(k,:));
            G0great0(k,:)=g0great0(k,:)+J.*Sz2.*(1i./(2*pi)).*(1-fermi(k,:)).*(E01(k,:)+F01(k,:)+E10(k,:)+F10(k,:));

            G1xless0(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sx2.*(E00(k,:)+F00(k,:))-1i.*Sy2.*(E10(k,:)+F10(k,:))+1i.*Sy2.*(E01(k,:)+F01(k,:))+1i.*Sx2.*(E11(k,:)+F11(k,:)));
            G1xgreat0(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sx2.*(E00(k,:)+F00(k,:))-1i.*Sy2.*(E10(k,:)+F10(k,:))+1i.*Sy2.*(E01(k,:)+F01(k,:))+1i.*Sx2.*(E11(k,:)+F11(k,:)));

            G1yless0(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sy2.*(E00(k,:)+F00(k,:))+1i.*Sx2.*(E10(k,:)+F10(k,:))-1i.*Sx2.*(E01(k,:)+F01(k,:))+1i.*Sy2.*(E11(k,:)+F11(k,:)));
            G1ygreat0(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sy2.*(E00(k,:)+F00(k,:))+1i.*Sx2.*(E10(k,:)+F10(k,:))-1i.*Sx2.*(E01(k,:)+F01(k,:))+1i.*Sy2.*(E11(k,:)+F11(k,:)));

            g1less0(k,:)=(1i./(2*pi)).*fermi(k,:).*(gS(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+g0(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));
            g1great0(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(gS(k).*(A0(k,:).*B0(k,:)+A1(k,:).*B1(k,:))+g0(k).*(A0(k,:).*B1(k,:)+A1(k,:).*B0(k,:)));

            G1zless0(k,:)=g1less0(k,:)-J.*Sz2.*(1i./(2*pi)).*fermi(k,:).*(E00(k,:)+F00(k,:));
            G1zgreat0(k,:)=g1great0(k,:)+J.*Sz2.*(1i./(2*pi)).*(1-fermi(k,:)).*(E00(k,:)+F00(k,:));

            if t(j) < t0
                if tau(i) < t0
                    K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                    K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                else
                    K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t0-tau(i)));
                    K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t0-tau(i)));
                end
            elseif t(j) < t1
                if tau(i) < t0
                    K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t0));
                    K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t0));
                elseif tau(i) < t1
                    K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-tau(i)));
                    K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-tau(i)));
                else
                    K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t1));
                    K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t(j)-t1));
                end
            else
                if tau(i) < t1
                    K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t1-tau(i)));
                    K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*eV(1).*(t1-tau(i)));
                else
                    K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                    K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                end
            end


%             %Function of (t,t')

            for m=1:2
                if t(j) < t0
                        A2(m,:) = -(1./2).*(exp(-1i.*w.*t(j))./(w-eps(m)+1i.*g(m)));
                elseif t(j) < t1
                        A2(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*g(m)).*(t(j)-t0)-1i.*w.*t0).*(1./(w-eps(m)+1i.*g(m))-1./(w+eV(k)-eps(m)+1i*g(m)))...
                        +exp(-1i.*eV(k).*(t(j)-t0)-1i.*w.*t(j))./(w+eV(k)-eps(m)+1i*g(m)));
                else
                        A2(m,:) = -(1./2).*(exp(-1i.*(eps(m)-1i.*g(m)).*(t(j))).*(1./(w-eps(m)+1i.*g(m)).*(...
                        exp(1i.*(eps(m)-1i.*g(m)-w).*t0)+exp(1i.*(eps(m)-1i.*g(m)-w).*t(j))-exp(1i.*(eps(m)-1i.*g(m)-w).*t1))...
                        +1./(w+eV(k)-eps(m)+1i*g(m)).*(exp(1i.*(eps(m)-1i.*g(m)-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*g(m)-w).*t0))));
                end

                if tau(i) < t0
                        B2(m,:) = -(1./2).*(exp(1i.*w.*tau(i))./(w-eps(m)-1i.*g(m)));
                elseif tau(i) < t1
                        B2(m,:) = -(1./2).*(exp(-1i.*(eps(m)+1i.*g(m)).*(t0-tau(i))+1i.*w.*t0).*(1./(w-eps(m)-1i.*g(m))-1./(w+eV(k)-eps(m)-1i*g(m)))...
                        +exp(1i.*eV(k).*(tau(i)-t0)+1i.*w.*tau(i))./(w+eV(k)-eps(m)-1i*g(m)));
                else
                        B2(m,:) = -(1./2).*(exp(1i.*(eps(m)+1i.*g(m)).*(tau(i))).*(1./(w-eps(m)-1i.*g(m)).*(...
                        exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*tau(i))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
                        +1./(w+eV(k)-eps(m)-1i*g(m)).*(exp(-1i.*(eps(m)+1i.*g(m)-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t0))));
                end

                for l=1:2

                    if t(j) < t0
                            C2(m,l,:)=(1./4).*exp(-1i.*w.*t(j))./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l)));
                    elseif t(j) < t1
                            C2(m,l,:)=(1./4).*(exp(-1i.*w.*t0-1i.*(eps(m)-1i.*g(m)).*(t(j)-t0))./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l)))...
                            +(exp(-1i.*w.*t(j)-1i.*eV(k).*(t(j)-t0))-exp(-1i.*w.*t0-1i.*(eps(m)-1i.*g(m)).*(t(j)-t0)))./((w+eV(k)-eps(m)+1i*g(m)).*(w+eV(k)-eps(l)+1i*g(l))));
                    else
                            C2(m,l,:) = (1./4).*(exp(-1i.*(eps(m)-1i.*g(m)).*(t(j))).*(1./((w-eps(m)+1i*g(m)).*(w-eps(l)+1i*g(l))).*(...
                            exp(1i.*(eps(m)-1i.*g(m)-w).*t0)+exp(1i.*(eps(m)-1i.*g(m)-w).*t(j))-exp(1i.*(eps(m)-1i.*g(m)-w).*t1))...
                            +1./((w+eV(k)-eps(m)+1i*g(m)).*(w+eV(k)-eps(l)+1i*g(l))).*(exp(1i.*(eps(m)-1i.*g(m)-w).*t1-1i.*eV(k).*(t1-t0))-exp(1i.*(eps(m)-1i.*g(m)-w).*t0))));
                    end

                    if tau(i) < t0
                            D2(m,l,:)=(1./4).*(exp(1i.*w.*tau(i))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))));
                    elseif tau(i) < t1
                            D2(m,l,:)=(1./4).*(exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-tau(i)))./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l)))...
                            +(exp(1i.*w.*tau(i)+1i.*eV(k).*(tau(i)-t0))-exp(1i.*w.*t0-1i.*(eps(m)+1i.*g(m)).*(t0-tau(i))))./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))));
                    else
                            D2(m,l,:) = (1./4).*(exp(1i.*(eps(m)+1i.*g(m)).*(tau(i))).*(1./((w-eps(m)-1i*g(m)).*(w-eps(l)-1i*g(l))).*(...
                            exp(-1i.*(eps(m)+1i.*g(m)-w).*t0)+exp(-1i.*(eps(m)+1i.*g(m)-w).*tau(i))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t1))...
                            +1./((w+eV(k)-eps(m)-1i*g(m)).*(w+eV(k)-eps(l)-1i*g(l))).*(exp(-1i.*(eps(m)+1i.*g(m)-w).*t1+1i.*eV(k).*(t1-t0))-exp(-1i.*(eps(m)+1i.*g(m)-w).*t0))));
                    end
                end
            end

            A02(k,:)=A2(1,:)+A2(2,:);
            A12(k,:)=A2(1,:)-A2(2,:);

            B02(k,:)=B2(1,:)+B2(2,:);
            B12(k,:)=B2(1,:)-B2(2,:);

            C002(k,:)=C2(1,1,:)+C2(1,2,:)+C2(2,1,:)+C2(2,2,:);
            C012(k,:)=C2(1,1,:)-C2(1,2,:)+C2(2,1,:)-C2(2,2,:);
            C102(k,:)=C2(1,1,:)+C2(1,2,:)-C2(2,1,:)-C2(2,2,:);
            C112(k,:)=C2(1,1,:)-C2(1,2,:)-C2(2,1,:)+C2(2,2,:);

            D002(k,:)=D2(1,1,:)+D2(1,2,:)+D2(2,1,:)+D2(2,2,:);
            D012(k,:)=D2(1,1,:)-D2(1,2,:)+D2(2,1,:)-D2(2,2,:);
            D102(k,:)=D2(1,1,:)+D2(1,2,:)-D2(2,1,:)-D2(2,2,:);
            D112(k,:)=D2(1,1,:)-D2(1,2,:)-D2(2,1,:)+D2(2,2,:);

            E002(k,:)=g0(k).*(C002(k,:).*B02(k,:)+C012(k,:).*B12(k,:))+gS(k).*(C002(k,:).*B12(k,:)+C012(k,:).*B02(k,:));
            E012(k,:)=gS(k).*(C002(k,:).*B02(k,:)+C012(k,:).*B12(k,:))+g0(k).*(C002(k,:).*B12(k,:)+C012(k,:).*B02(k,:));
            E102(k,:)=g0(k).*(C102(k,:).*B02(k,:)+C112(k,:).*B12(k,:))+gS(k).*(C102(k,:).*B12(k,:)+C112(k,:).*B02(k,:));
            E112(k,:)=gS(k).*(C102(k,:).*B02(k,:)+C112(k,:).*B12(k,:))+g0(k).*(C102(k,:).*B12(k,:)+C112(k,:).*B02(k,:));
            F002(k,:)=g0(k).*(A02(k,:).*D002(k,:)+A12(k,:).*D012(k,:))+gS(k).*(A02(k,:).*D012(k,:)+A12(k,:).*D002(k,:));
            F012(k,:)=g0(k).*(A02(k,:).*D102(k,:)+A12(k,:).*D112(k,:))+gS(k).*(A02(k,:).*D112(k,:)+A12(k,:).*D102(k,:));
            F102(k,:)=gS(k).*(A02(k,:).*D002(k,:)+A12(k,:).*D012(k,:))+g0(k).*(A02(k,:).*D012(k,:)+A12(k,:).*D002(k,:));
            F112(k,:)=gS(k).*(A02(k,:).*D102(k,:)+A12(k,:).*D112(k,:))+g0(k).*(A02(k,:).*D112(k,:)+A12(k,:).*D102(k,:));

            g0less02(k,:)=(1i./(2*pi)).*fermi(k,:).*(g0(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+gS(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));
            g0great02(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(g0(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+gS(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));

            G0less02(k,:)=g0less02(k,:)-J.*Sz2.*(1i./(2*pi)).*fermi(k,:).*(E012(k,:)+F012(k,:)+E102(k,:)+F102(k,:));
            G0great02(k,:)=g0great02(k,:)+J.*Sz2.*(1i./(2*pi)).*(1-fermi(k,:)).*(E012(k,:)+F012(k,:)+E102(k,:)+F102(k,:));

            G1xless02(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sx2.*(E002(k,:)+F002(k,:))-1i.*Sy2.*(E102(k,:)+F102(k,:))+1i.*Sy2.*(E012(k,:)+F012(k,:))+1i.*Sx2.*(E112(k,:)+F112(k,:)));
            G1xgreat02(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sx2.*(E002(k,:)+F002(k,:))-1i.*Sy2.*(E102(k,:)+F102(k,:))+1i.*Sy2.*(E012(k,:)+F012(k,:))+1i.*Sx2.*(E112(k,:)+F112(k,:)));

            G1yless02(k,:)=-J.*(1i./(2*pi)).*fermi(k,:).*(Sy2.*(E002(k,:)+F002(k,:))+1i.*Sx2.*(E102(k,:)+F102(k,:))-1i.*Sx2.*(E012(k,:)+F012(k,:))+1i.*Sy2.*(E112(k,:)+F112(k,:)));
            G1ygreat02(k,:)=J.*(1i./(2*pi)).*(1-fermi(k,:)).*(Sy2.*(E002(k,:)+F002(k,:))+1i.*Sx2.*(E102(k,:)+F102(k,:))-1i.*Sx2.*(E012(k,:)+F012(k,:))+1i.*Sy2.*(E112(k,:)+F112(k,:)));

            g1less02(k,:)=(1i./(2*pi)).*fermi(k,:).*(gS(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+g0(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));
            g1great02(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*(gS(k).*(A02(k,:).*B02(k,:)+A12(k,:).*B12(k,:))+g0(k).*(A02(k,:).*B12(k,:)+A12(k,:).*B02(k,:)));

            G1zless02(k,:)=g1less02(k,:)-J.*Sz2.*(1i./(2*pi)).*fermi(k,:).*(E002(k,:)+F002(k,:));
            G1zgreat02(k,:)=g1great02(k,:)+J.*Sz2.*(1i./(2*pi)).*(1-fermi(k,:)).*(E002(k,:)+F002(k,:));

            %No change in voltage and no spin dependency
            %greensfunctionAlt1
            %No change in voltage
            %greensfunctionAlt2
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

    G0greattotAlt(j,:)=G0great;
    G1xgreattotAlt(j,:)=G1xgreat;
    G1ygreattotAlt(j,:)=G1ygreat;
    G1zgreattotAlt(j,:)=G1zgreat;

    jH=1i.*J.^2.*(G0less.*G0great2-G0great.*G0less2-G1xless.*G1xgreat2+G1xgreat.*G1xless2...
        -G1yless.*G1ygreat2+G1ygreat.*G1yless2-G1zless.*G1zgreat2+G1zgreat.*G1zless2);
    jDMx=-J.^2.*(G0less.*G1xgreat2-G0great.*G1xless2-G1xless.*G0great2+G1xgreat.*G0less2);
    jDMy=-J.^2.*(G0less.*G1ygreat2-G0great.*G1yless2-G1yless.*G0great2+G1ygreat.*G0less2);
    jDMz=-J.^2.*(G0less.*G1zgreat2-G0great.*G1zless2-G1zless.*G0great2+G1zgreat.*G0less2);
    jIxx=2*1i.*J.^2.*(G1xless.*G1xgreat2-G1xgreat.*G1xless2);
    jIyy=2*1i.*J.^2.*(G1yless.*G1ygreat2-G1ygreat.*G1yless2);
    jIzz=2*1i.*J.^2.*(G1zless.*G1zgreat2-G1zgreat.*G1zless2);
    jIxy=1i.*J.^2.*(G1xless.*G1ygreat2-G1xgreat.*G1yless2+G1yless.*G1xgreat2-G1ygreat.*G1xless2);
    jIyx=jIxy;
    jIxz=1i.*J.^2.*(G1xless.*G1zgreat2-G1xgreat.*G1zless2+G1zless.*G1xgreat2-G1zgreat.*G1xless2);
    jIzx=jIxz;
    jIyz=1i.*J.^2.*(G1yless.*G1zgreat2-G1ygreat.*G1zless2+G1zless.*G1ygreat2-G1zgreat.*G1yless2);
    jIzy=jIyz;

    ejxbare=2*1i*J*((epsilon*G0less+wL/2.*G1zless).*G1xgreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1xless2...
        +(epsilon*G1xless-1i.*G1yless).*G0great2-(epsilon*G1xgreat-1i.*G1ygreat).*G0less2...
        -1i.*(epsilon*G1yless+1i.*G1xless).*G1zgreat2+1i.*(epsilon*G1zless+wL/2).*G1ygreat2...
        +1i.*(epsilon*G1ygreat+1i.*G1xgreat).*G1zless2-1i.*(epsilon*G1zgreat+wL/2).*G1yless2);
    ejybare=2*1i*J*((epsilon*G0less+wL/2.*G1zless).*G1ygreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1yless2...
        +(epsilon*G1yless+1i.*G1xless).*G0great2-(epsilon*G1ygreat+1i.*G1xgreat).*G0less2...
        -1i.*(epsilon*G1zless+wL/2).*G1xgreat2+1i.*(epsilon*G1xless-1i.*G1yless).*G1zgreat2...
        +1i.*(epsilon*G1zgreat+wL/2).*G1xless2-1i.*(epsilon*G1xgreat-1i.*G1ygreat).*G1zless2);
    ejzbare=2*1i*J*((epsilon*G0less+wL/2.*G1zless).*G1zgreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1zless2...
        +(epsilon*G1zless+wL/2).*G0great2-(epsilon*G1zgreat2+wL/2).*G0less2...
        -1i.*(epsilon*G1xless-1i.*G1yless).*G1ygreat2+1i.*(epsilon*G1yless+1i.*G1xless).*G1xgreat2...
        +1i.*(epsilon*G1xgreat-1i.*G1ygreat).*G1yless2-1i.*(epsilon*G1ygreat+1i.*G1xgreat).*G1xless2);

    ejx=trapz(tau,ejxbare);
    ejy=trapz(tau,ejybare);
    ejz=trapz(tau,ejzbare);

    %�ndrat s� att den ber�knar den nya ekvationen med Sx2, Sy2 och Sz2 och SxT2, SyT2, SzT2
    dSxbare=jH.*(Sy2.*SzT2-Sz2.*SyT2)+Sy2.*(jIzx.*SxT2+jIzy.*SyT2+jIzz.*SzT2)-Sz2.*(jIyx.*SxT2+jIyy.*SyT2+jIyz.*SzT2)-Sy2.*(jDMx.*SyT2-jDMy.*SxT2)+Sz2.*(jDMz.*SxT2-jDMx.*SzT2);
    dSybare=jH.*(Sz2.*SxT2-Sx2.*SzT2)+Sz2.*(jIxx.*SxT2+jIxy.*SyT2+jIxz.*SzT2)-Sx2.*(jIzx.*SxT2+jIzy.*SyT2+jIzz.*SzT2)-Sz2.*(jDMy.*SzT2-jDMz.*SyT2)+Sx2.*(jDMx.*SyT2-jDMy.*SxT2);
    dSzbare=jH.*(Sx2.*SyT2-Sy2.*SxT2)+Sx2.*(jIyx.*SxT2+jIyy.*SyT2+jIyz.*SzT2)-Sy2.*(jIxx.*SxT2+jIxy.*SyT2+jIxz.*SzT2)-Sx2.*(jDMz.*SxT2-jDMx.*SzT2)+Sy2.*(jDMy.*SzT2-jDMz.*SyT2);
    dSx2=-(wL+ejz).*Sy2+ejy.*Sz2+trapz(tau,dSxbare);
    dSy2=-ejx.*Sz2+(wL+ejz).*Sx2+trapz(tau,dSybare);
    dSz2=-ejy.*Sx2+ejx.*Sy2+trapz(tau,dSzbare);

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
        SxNonorm(j+1)=Sx(j)+tstep*real(dSx);
        SyNonorm(j+1)=Sy(j)+tstep*real(dSy);
        SzNonorm(j+1)=Sz(j)+tstep*real(dSz);

        Svector=[SxNonorm(j+1),SyNonorm(j+1),SzNonorm(j+1)];
        Slength(j)=norm(Svector);
        SvectorNorm=Svector/norm(Svector);
        Sx(j+1)=SvectorNorm(1);
        Sy(j+1)=SvectorNorm(2);
        Sz(j+1)=SvectorNorm(3);
    end
end


%In SI units
%tconv=6.5821189*10^(-16)/10^-3;%in s, hbar in eVs/meV
%t=tconv.*t;
Iconv=1.602176565*10^(-19)/(6.58211928*10^(-16))*10^-3;%in A, elementary charge divided by hbar in eVs times 10^-3

Ic=(Ic0+Ic1).*Iconv;
Isx=(Isx0+Isx1).*Iconv;
Isy=(Isy0+Isy1).*Iconv;
Isz=(Isz0+Isz1).*Iconv;
