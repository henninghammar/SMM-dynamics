% Polarized single-molecule magnet with pulse
% ================================================
% Author: Henning Hammar
% ------------------
% Calculation of spin dynamics and current through a polarized single-molecule magnet
% for a pulse.
%
% This file set all values needed in order to perform a run

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
J0=0.1*gamma; %Coupling strength
Sz0=cos(pi/4); %Spin z-component
Sxy=sin(pi/4); %Spin xy-component
wL=gfactor*myB*B; %Frequency
epsilon=0; %Energy level of the quantum dot
eps(1)=epsilon+0.5*wL;
eps(2)=epsilon-0.5*wL;

%Time variables and time and energy step-size
tmax=10;
tstep=0.1;
tstep2=0.1/gamma;
tback=200*gamma;
t0=0;%50*tstep;
t1=3;

kB=8.617324*10^-2; %Boltzmanns constant, in meV*K^-1
T(1)=1; %Temperature in K
T(2)=1;

main

toc
%
% clear G0less G0great G1xless G1xgreat G1yless G1ygreat G1zless G1zgreat G0less2 G0great2 G1xless2 G1xgreat2 G1yless2 G1ygreat2 G1zless2 G1zgreat2 Kless Kgreat
% clear G0greattot   G1xgreattot   G1ygreattot  G1zgreattot   G0greattot2   G1xgreattot2    G1ygreattot2    G1zgreattot2
% clear G0greattotAlt    G1xgreattotAlt    G1ygreattotAlt    G1zgreattotAlt
%
% clear dSxbarejH dSybarejH dSzbarejH dSxbarejDMx dSybarejDMy dSzbarejDMz dSxbarejIx dSybarejIy dSzbarejIz
% clear jHtot jDMxtot jDMytot jDMztot jIxxtot jIyytot jIzztot jIxytot jIxztot jIyztot
% clear    barejDMfieldx    barejDMfieldy    barejDMfieldz    barejHfield barejIfieldx barejIfieldy barejIfieldz

%Savedata
outputFolder = 'output';
outputFilename = sprintf('%s/test1.mat', outputFolder);
save(outputFilename)
