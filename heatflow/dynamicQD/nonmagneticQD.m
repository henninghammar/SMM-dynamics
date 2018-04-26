%matlabpool open

clear

addpath('../../src')

tic

%Initial values
kB=8.6173324*10^-2; %Boltzmanns constant, in meV*K^-1
T=[1 1]; %Temperature in K
beta=1/(kB*T(1)); %Beta value
delta=10^-6; %i*delta, size of delta
g0L=1;
g0R=1;
g0=1; %g0 in eV (Coupling strength between leads and the dot)
eVL=1; %bias voltage on left lead
eVR=-1; %bias voltage on right lead
mu = [eVL; eVR];
eps=0; %Energy level of the quantum dot

%Time variables and time and energy step-size
t0=0;
t1=100;
tmax=1;
tstep=0.1;
tstep2=0.01;
tback=500;
t=[t0:tstep:tmax];
step=0.1;
w=[-50:step:50];

%Run the program
fermi=1./(1+exp(beta.*w));
fermi=[fermi; fermi];
Ic=zeros(1,length(t));

for j=1:length(t)
    tau = [-tstep2*tback+t(j):tstep2:t(j)];

    timestep=t(j);
    timeused=toc;
    disp(['Timestep: ' num2str(timestep)])
    disp(['Timeused: ' num2str(timeused)])

    for i=1:tback+1

      %Green's function for a voltage pulse
      for k=1:2
          if tau(i) < t0
                  A = -(1./2).*(exp(-1i.*w.*tau(i))./(w-eps+1i.*g0/2));
          elseif tau(i) < t1
                  A = -(1./2).*(exp(-1i.*(eps-1i.*g0/2).*(tau(i)-t0)-1i.*w.*t0).*(1./(w-eps+1i.*g0/2)-1./(w+mu(k)-eps+1i*g0/2))...
                  +exp(-1i.*mu(k).*(tau(i)-t0)-1i.*w.*tau(i))./(w+mu(k)-eps+1i*g0/2));
          else
                  A = -(1./2).*(exp(-1i.*(eps-1i.*g0/2).*(tau(i))).*(1./(w-eps+1i.*g0/2).*(...
                  exp(1i.*(eps-1i.*g0/2-w).*t0)+exp(1i.*(eps-1i.*g0/2-w).*tau(i))-exp(1i.*(eps-1i.*g0/2-w).*t1))...
                  +1./(w+mu(k)-eps+1i*g0/2).*(exp(1i.*(eps-1i.*g0/2-w).*t1-1i.*mu(k).*(t1-t0))-exp(1i.*(eps-1i.*g0/2-w).*t0))));
          end

          if t(j) < t0
                  B = -(1./2).*(exp(1i.*w.*t(j))./(w-eps-1i.*g0/2));
          elseif t(j) < t1
                  B = -(1./2).*(exp(-1i.*(eps+1i.*g0/2).*(t0-t(j))+1i.*w.*t0).*(1./(w-eps-1i.*g0/2)-1./(w+mu(k)-eps-1i*g0/2))...
                  +exp(1i.*mu(k).*(t(j)-t0)+1i.*w.*t(j))./(w+mu(k)-eps-1i*g0/2));
          else
                  B = -(1./2).*(exp(1i.*(eps+1i.*g0/2).*(t(j))).*(1./(w-eps-1i.*g0/2).*(...
                  exp(-1i.*(eps+1i.*g0/2-w).*t0)+exp(-1i.*(eps+1i.*g0/2-w).*t(j))-exp(-1i.*(eps+1i.*g0/2-w).*t1))...
                  +1./(w+mu(k)-eps-1i*g0/2).*(exp(-1i.*(eps+1i.*g0/2-w).*t1+1i.*mu(k).*(t1-t0))-exp(-1i.*(eps+1i.*g0/2-w).*t0))));
          end

          G0less0(k,:)=(1i./(2*pi)).*fermi(k,:).*g0.*A.*B;
          G0great0(k,:)=-(1i./(2*pi)).*(1-fermi(k,:)).*g0.*A.*B;
      end

      %Integrating over energies for the different Green's functions

      G0less(i)=trapz(w,G0less0(1,:)+G0less0(2,:));
      G0great(i)=trapz(w,G0great0(1,:)+G0great0(2,:));

      selfenergyK
      %[energyKLless(i), energyKLgreat(i), energyKRless(i), energyKRgreat(i)] = energyselfenergyK(w2, t(j), tau(i), t0, t1, mu, T, kB);

      for k=1:2
          if t(j) < t0
              if tau(i) < t0
                  energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                  energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
              else
                  energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t0-tau(i)));
                  energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t0-tau(i)));
              end
          elseif t(j) < t1
              if tau < t0
                  energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))-1i.*mu(k).*(t(j)-t0));
                  energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))-1i.*mu(k).*(t(j)-t0));
              elseif tau < t1
                  energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))-1i.*mu(k).*(t(j)-tau(i)));
                  energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))-1i.*mu(k).*(t(j)-tau(i)));
              else
                  energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))-1i.*mu(k).*(t(j)-t1));
                  energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))-1i.*mu(k).*(t(j)-t1));
              end
          else
              if tau < t1
                  energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t1-tau(i)));
                  energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i))).*exp(-1i.*mu(k).*(t1-tau(i)));
              else
                  energyK0less(k,:)=w.*fermi(k,:)./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
                  energyK0great(k,:)=w.*(1-fermi(k,:))./(2*pi).*exp(-1i.*w.*(t(j)-tau(i)));
              end
          end
      end

      energyKLless(i)=trapz(w,energyK0less(1,:));
      energyKLgreat(i)=trapz(w,energyK0great(1,:));
      energyKRless(i)=trapz(w,energyK0less(2,:));
      energyKRgreat(i)=trapz(w,energyK0great(2,:));

      energycurrentLbare(i)=trapz(w,(G0less0(1,:)+G0less0(2,:)).*energyK0great(1,:)+(G0great0(1,:)+G0great0(2,:)).*energyK0less(1,:));
      energycurrentRbare(i)=trapz(w,(G0less0(1,:)+G0less0(2,:)).*energyK0great(2,:)+(G0great0(1,:)+G0great0(2,:)).*energyK0less(2,:));

    end

     IcL0=G0less.*KLgreat+G0great.*KLless;
     IcL(j)=-4*g0L*imag(trapz(tau,IcL0));
     IcR0=G0less.*KRgreat+G0great.*KRless;
     IcR(j)=-4*g0R*imag(trapz(tau,IcR0));
     %
     IeLbare=G0less.*energyKLgreat+G0great.*energyKLless;
     IeL(j)=-4.*g0L*imag(trapz(tau,IeLbare));
     IeRbare=G0less.*energyKRgreat+G0great.*energyKRless;
     IeR(j)=-4.*g0R*imag(trapz(tau,IeRbare));
     %IeL(j)=-4.*g0L*imag(trapz(tau,energycurrentLbare));
     %IeR(j)=-4.*g0R*imag(trapz(tau,energycurrentRbare));

     IqL(j)=IeL(j)-mu(1)*IcL(j);
     IqR(j)=IeR(j)-mu(2)*IcR(j);
end

toc

outputFolder = 'output';
outputFilename = sprintf('%s/test23.mat', outputFolder);
save(outputFilename)

figure
plot(t, IcL)

figure
plot(t, IeL)

figure
plot(t, IqL)
