%Green's function for a voltage pulse with t,t' ordering

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