function [energyKLless, energyKLgreat, energyKRless, energyKRgreat] = energyselfenergyK(w, t, tau, t0, t1, mu, T, kB)
%Calculate the K integral of the self-energy for different times considering
%a voltage pulse

fermifunction

for k=1:2
    if t < t0
        if tau < t0
            energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t-tau));
            energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t-tau));
        else
            energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t0-tau));
            energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t0-tau));
        end
    elseif t < t1
        if tau < t0
            energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t-t0));
            energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t-t0));
        elseif tau < t1
            energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t-tau));
            energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t-tau));
        else
            energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t-t1));
            energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t-t1));
        end
    else
        if tau < t1
            energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t1-tau));
            energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t-tau)).*exp(-1i.*mu(k).*(t1-tau));
        else
            energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t-tau));
            energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t-tau));
        end
    end
end

energyKLless=trapz(w,energyK0less(1,:));
energyKLgreat=trapz(w,energyK0great(1,:));
energyKRless=trapz(w,energyK0less(2,:));
energyKRgreat=trapz(w,energyK0great(2,:));

end
