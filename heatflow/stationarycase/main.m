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
eV = gamma; %bias voltage
mu(1)=eV/2; %chemical potential on left lead
mu(2)=-eV/2; %chemical potential on right lead
J=1;%0.5*gamma; %Coupling strength

gfactor=2; %g-factor
myB=5.78838175*10^(-2); %Bohr magneton in meV*T^-1
B=1; %Magnetic field in Tesla
wL=gfactor*myB*B; %Frequency
epsilon=0.5; %Energy level of the quantum dot
eps(1)=epsilon+0.5*wL;
eps(2)=epsilon-0.5*wL;

S = [0,0,1];

kB=8.617324*10^-2; %Boltzmanns constant, in meV*K^-1
T(1)=1; %Temperature in K
T(2)=1;
beta(1)=1/(kB*T(1));
beta(2)=1/(kB*T(2));

dw = 0.1;
w=[-20:dw:20];

eVvector = [-10:0.1:10]./gamma;
for i = 1:length(eVvector)
  munew(1) = eVvector(i)/2;
  munew(2) = -eVvector(i)/2;
  [Icbias(i), Isxbias(i), Isybias(i), Iszbias(i), Iebias(i), Iqbias(i), DOSbias(i,:), MDOSzbias(i,:)] = stationarycurrent(pL, pR, gamma, munew, eps, epsilon, w, dw, J, S, wL, beta);
  %[JH, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Dx, Dy, Dz, ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, beta);
end

epsvector = [-10:0.1:10]./gamma;
for i = 1:length(epsvector)
  epsnew(1)=epsvector(i)+0.5*wL;
  epsnew(2)=epsvector(i)-0.5*wL;
  [Icgate(i), Isxgate(i), Isygate(i), Iszgate(i), Iegate(i), Iqgate(i), DOSgate(i,:), MDOSzgate(i,:)] = stationarycurrent(pL, pR, gamma, mu, epsnew, epsilon, w, dw, J, S, wL, beta);
  %[JH, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Dx, Dy, Dz, ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, beta);
end

Jvector = [0:0.1:10]./gamma;
for i = 1:length(Jvector)
  [IcJ(i), IsxJ(i), IsyJ(i), IszJ(i), IeJ(i), IqJ(i), DOSJ(i,:), MDOSzJ(i,:)] = stationarycurrent(pL, pR, gamma, mu, eps, epsilon, w, dw, Jvector(i), S, wL, beta);
  %[JH, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Dx, Dy, Dz, ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, beta);
end

tempvector = [1:0.1:10];
for i = 1:length(tempvector)
  Tnew(1)=1; %Temperature in K
  Tnew(2)=tempvector(i);
  betanew(1)=1/(kB*Tnew(1));
  betanew(2)=1/(kB*Tnew(2));
  [Ictemp(i), Isxtemp(i), Isytemp(i), Isztemp(i), Ietemp(i), Iqtemp(i), DOStemp(i,:), MDOSztemp(i,:)] = stationarycurrent(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, betanew);
  %[JH, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Dx, Dy, Dz, ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, beta);
end

toc

plotstationary

%Savedata
outputFolder = 'output';
outputFilename = sprintf('%s/test5.mat', outputFolder);
save(outputFilename)
