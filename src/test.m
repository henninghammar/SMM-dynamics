%Testing different Green's functions

%Initial values
pL=0; %Polarization of gamma_up and gamma_down
pR=0;
gamma=1;
g0(1)=gamma;
g0(2)=gamma;
eV(1)=1; %bias voltage on left lead
eV(2)=-1; %bias voltage on right lead
gfactor=2; %g-factor
myB=5.78838175*10^(-2); %Bohr magneton in meV*T^-1
B=1; %Magnetic field in Tesla
J=0.1; %Coupling strength
Sz0=cos(pi/4); %Spin z-component
Sxy=sin(pi/4); %Spin xy-component
wL=gfactor*myB*B; %Frequency
epsilon=0; %Energy level of the quantum dot
eps(1)=epsilon;%+0.5*wL;
eps(2)=epsilon;%-0.5*wL;

tau = 2;
t = 4;
t0=0;%50*tstep;
t1=3;

kB=8.617324*10^-2; %Boltzmanns constant, in meV*K^-1
T(1)=1; %Temperature in K
T(2)=1;

step=0.1;
w=[-20:step:20];

%Defining a vector for the spin
Sx=-Sxy*sin(wL*t0);
Sy=Sxy*cos(wL*t0);
Sz=Sz0;
S = [Sx, Sy, Sz];

%Polarized gamma for leads with different magnetization
polarizedgamma

%Defining fermi-functions
fermifunction

[G0less, G0great, G1xless, G1xgreat, G1yless, G1ygreat, G1zless, G1zgreat] = greensfunction(tau, t, t0, t1, eps, g, g0, gS, eV, w, fermi, S, J);
[G0less2, G0great2, G1xless2, G1xgreat2, G1yless2, G1ygreat2, G1zless2, G1zgreat2] = SMMnonpolarizedgreensfunction(tau, t, t0, t1, epsilon, g0tot, g0, eV, w, fermi, S, J);
