% Test of the greens function
% ================================================
% Author: Henning Hammar
% ------------------

clear

addpath('../../src')

tic

%Initial values
pL=0; %Polarization of gamma_up and gamma_down
pR=0;
gamma=1;
g0(1)=gamma;
g0(2)=gamma;
eV(1)=gamma; %bias voltage on left lead
eV(2)=-gamma; %bias voltage on right lead
gfactor=2; %g-factor
myB=5.78838175*10^(-2); %Bohr magneton in meV*T^-1
B=1; %Magnetic field in Tesla
J=0.1*gamma; %Coupling strength
S = [0, sin(pi/4), cos(pi/4)];
wL=gfactor*myB*B; %Frequency
epsilon=0; %Energy level of the quantum dot
eps(1)=epsilon+0.5*wL;
eps(2)=epsilon-0.5*wL;

%Time variables and time and energy step-size
tstep=0.1;
tstep2=0.1/gamma;
tback=200*gamma;
t0=0;%50*tstep;
t1=3;

kB=8.617324*10^-2; %Boltzmanns constant, in meV*K^-1
T(1)=1; %Temperature in K
T(2)=1;

%Integration values in time and energy
t=5;
step=0.1;
w=[-20:step:20];

%Polarized gamma for leads with different magnetization
polarizedgamma

%Defining fermi-functions
fermifunction

%Initiating some variables for speed
%initiatingvariables

tau = [-tstep2*tback+t:tstep2:t];

for i=1:tback+1
    %Green's function of (t,tau) for each timestep tau with integration over energies w/omega
    [G0less(i), G0great(i), G1xless(i), G1xgreat(i), G1yless(i), G1ygreat(i), G1zless(i), G1zgreat(i)] = greensfunction(t, tau(i), t0, t1, eps, g, g0, gS, eV, w, fermi, S, J);
end

%Savedata
outputFolder = 'output';
outputFilename = sprintf('%s/test6.mat', outputFolder);
%save(outputFilename)
