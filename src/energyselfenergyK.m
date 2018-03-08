%Calculate the K integral of the self-energy for different times considering
%a voltage pulse

for m=1:2
  for k=1:2
      if t(j) < t0
          if tau(i) < t0
              energyK0less(m,k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
              energyK0great(m,k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
          else
              energyK0less(m,k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t0-tau(i)));
              energyK0great(m,k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t0-tau(i)));
          end
      elseif t(j) < t1
          if tau(i) < t0
              energyK0less(m,k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t(j)-t0));
              energyK0great(m,k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t(j)-t0));
          elseif tau(i) < t1
              energyK0less(m,k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t(j)-tau(i)));
              energyK0great(m,k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t(j)-tau(i)));
          else
              energyK0less(m,k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t(j)-t1));
              energyK0great(m,k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t(j)-t1));
          end
      else
          if tau(i) < t1
              energyK0less(m,k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t1-tau(i)));
              energyK0great(m,k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(m).*(t1-tau(i)));
          else
              energyK0less(m,k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
              energyK0great(m,k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
          end
      end
  end
end

energyKLless(i)=trapz(w,energyK0less(1,1,:)+energyK0less(1,2,:));
energyKLgreat(i)=trapz(w,energyK0great(1,1,:)+energyK0great(1,2,:));
energyKRless(i)=trapz(w,energyK0less(2,1,:)+energyK0less(2,2,:));
energyKRgreat(i)=trapz(w,energyK0great(2,1,:)+energyK0great(2,2,:));
