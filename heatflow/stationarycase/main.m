% Stationary currents
% ================================================
% Author: Henning Hammar
% ------------------
% Calculation of the current and other observables for the stationary case of an SMM in a tunnel junction.

clear

addpath('../../src')

tic

%Initial values
pL=0; %Polarization of gamma_up and gamma_down
pR=0;
gamma=1;
eV = 0;%gamma; %bias voltage
mu(1)=eV/2; %chemical potential on left lead
mu(2)=-eV/2; %chemical potential on right lead
J=1;%0.5*gamma; %Coupling strength

gfactor=2; %g-factor
myB=5.78838175*10^(-2); %Bohr magneton in meV*T^-1
B=0; %Magnetic field in Tesla
wL=gfactor*myB*B; %Frequency
epsilon=0; %Energy level of the quantum dot
eps(1)=epsilon+0.5*wL;
eps(2)=epsilon-0.5*wL;

S = [0,0,1];

kB=8.617324*10^-2; %Boltzmanns constant, in meV*K^-1
T(1)=1; %Temperature in K
T(2)=1;
beta(1)=1/(kB*T(1));
beta(2)=1/(kB*T(2));

dw = 0.1;
w=[-10:dw:10];

eVvector = [-10:0.1:10]./gamma;
for i = 1:length(eVvector)
  munew(1) = eVvector(i)/2;
  munew(2) = -eVvector(i)/2;
  [Icbias(i), Isxbias(i), Isybias(i), Iszbias(i), Iebias(i), Iqbias(i), DOSbias(i,:), MDOSzbias(i,:)] = stationarycurrent(pL, pR, gamma, munew, eps, epsilon, w, dw, J, S, wL, beta);
  [JHbias(i), Isingbias(1,i), Isingbias(2,i), Isingbias(3,i), Isingbias(4,i), Isingbias(5,i), Isingbias(6,i), DMbias(1,i), DMbias(2,i), DMbias(3,i), ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, munew, eps, epsilon, w, dw, J, S, wL, beta);
end

epsvector = [-10:0.1:10]./gamma;
for i = 1:length(epsvector)
  epsnew(1)=epsvector(i)+0.5*wL;
  epsnew(2)=epsvector(i)-0.5*wL;
  [Icgate(i), Isxgate(i), Isygate(i), Iszgate(i), Iegate(i), Iqgate(i), DOSgate(i,:), MDOSzgate(i,:)] = stationarycurrent(pL, pR, gamma, mu, epsnew, epsilon, w, dw, J, S, wL, beta);
  [JHgate(i), Isinggate(1,i), Isinggate(2,i), Isinggate(3,i), Isinggate(4,i), Isinggate(5,i), Isinggate(6,i), DMgate(1,i), DMgate(2,i), DMgate(3,i), ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, epsnew, epsilon, w, dw, J, S, wL, beta);
end

Jvector = [0:0.1:10]./gamma;
for i = 1:length(Jvector)
  [IcJ(i), IsxJ(i), IsyJ(i), IszJ(i), IeJ(i), IqJ(i), DOSJ(i,:), MDOSzJ(i,:)] = stationarycurrent(pL, pR, gamma, mu, eps, epsilon, w, dw, Jvector(i), S, wL, beta);
  [JHJ(i), IsingJ(1,i), IsingJ(2,i), IsingJ(3,i), IsingJ(4,i), IsingJ(5,i), IsingJ(6,i), DMJ(1,i), DMJ(2,i), DMJ(3,i), ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, eps, epsilon, w, dw, Jvector(i), S, wL, beta);
end

tempvector = [1:0.1:10];
for i = 1:length(tempvector)
  Tnew(1)=1; %Temperature in K
  Tnew(2)=tempvector(i);
  betanew(1)=1/(kB*Tnew(1));
  betanew(2)=1/(kB*Tnew(2));
  [Ictemp(i), Isxtemp(i), Isytemp(i), Isztemp(i), Ietemp(i), Iqtemp(i), DOStemp(i,:), MDOSztemp(i,:)] = stationarycurrent(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, betanew);
  [JHtemp(i), Isingtemp(1,i), Isingtemp(2,i), Isingtemp(3,i), Isingtemp(4,i), Isingtemp(5,i), Isingtemp(6,i), DMtemp(1,i), DMtemp(2,i), DMtemp(3,i), ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, betanew);
end

toc

%Savedata
outputFolder = 'output';
outputFilename = sprintf('%s/paranobianogateSz1.mat', outputFolder);
save(outputFilename)

plotstationary
