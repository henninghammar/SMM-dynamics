%This file set all values needed in order to perform a run
%It can be altered for different projects and different systems

clear
tic

%Initial values
pL=0; %Polarization of gamma_up and gamma_down
pR=0;
gamma=1;
g0(1)=gamma;
g0(2)=gamma;
eV(1)=4; %bias voltage on left lead
eV(2)=-4; %bias voltage on right lead
gfactor=2; %g-factor
myB=5.78838175*10^(-2); %Bohr magneton in meV*T^-1
B=0; %Magnetic field in Tesla
J0=1.2; %Coupling strength
Sz0=cos(pi/4); %Spin z-component
Sxy=sin(pi/4); %Spin xy-component
wL=gfactor*myB*B; %Frequency
epsilon=0; %Energy level of the quantum dot
eps(1)=epsilon+0.5*wL;
eps(2)=epsilon-0.5*wL;

%Time variables and time and energy step-size
tscale=1;
tmax=0.2;
tstep=0.05;
tback=200;
t0=0;%50*tstep;
t1=3;

kB=8.617324*10^-2; %Boltzmanns constant, in meV*K^-1
T(1)=1; %Temperature in K
T(2)=1;

main

toc

clear G0less G0great G1xless G1xgreat G1yless G1ygreat G1zless G1zgreat G0less2 G0great2 G1xless2 G1xgreat2 G1yless2 G1ygreat2 G1zless2 G1zgreat2 Kless Kgreat   
clear G0greattot   G1xgreattot   G1ygreattot  G1zgreattot   G0greattot2   G1xgreattot2    G1ygreattot2    G1zgreattot2
clear G0greattotAlt    G1xgreattotAlt    G1ygreattotAlt    G1zgreattotAlt

clear dSxbarejH dSybarejH dSzbarejH dSxbarejDMx dSybarejDMy dSzbarejDMz dSxbarejIx dSybarejIy dSzbarejIz
clear jHtot jDMxtot jDMytot jDMztot jIxxtot jIyytot jIzztot jIxytot jIxztot jIyztot
clear    barejDMfieldx    barejDMfieldy    barejDMfieldz    barejHfield barejIfieldx barejIfieldy barejIfieldz

% Savedata
%outputFolder = 'outfolder';
%outputFilename = sprintf('%s/pulsebiasnha%02d.mat', outputFolder, bias);
%save(outputFilename)