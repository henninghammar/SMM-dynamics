%Calculate the K integral of the sel-energy for different times considering
%a voltage pulse

for k=1:2
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
end

Kless(i)=trapz(w,K0less(1,:)+K0less(2,:));
Kgreat(i)=trapz(w,K0great(1,:)+K0great(2,:));