%Calculate the K integral of the self-energy for different times considering
%a voltage pulse

for k=1:2
    if t(j) < t0
        if tau(i) < t0
            K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
            K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
        else
            K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t0-tau(i)));
            K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t0-tau(i)));
        end
    elseif t(j) < t1
        if tau(i) < t0
            K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t(j)-t0));
            K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t(j)-t0));
        elseif tau(i) < t1
            K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t(j)-tau(i)));
            K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t(j)-tau(i)));
        else
            K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t(j)-t1));
            K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t(j)-t1));
        end
    else
        if tau(i) < t1
            K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t1-tau(i)));
            K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t1-tau(i)));
        else
            K0less(k,:)=fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
            K0great(k,:)=(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
        end
    end
end

KLless(i)=trapz(w,K0less(1,:));
KLgreat(i)=trapz(w,K0great(1,:));
KRless(i)=trapz(w,K0less(2,:));
KRgreat(i)=trapz(w,K0great(2,:));
