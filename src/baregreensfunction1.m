%Green's function for a voltage pulse with t,t' ordering

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
            