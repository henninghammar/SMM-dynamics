% Polarized single-molecule magnet with pulse with tdLLG with inertia
% ================================================
% Author: Henning Hammar
% ------------------
% Calculation of the spin dynamics and current through a polarized single-molecule magnet
% for a pulse using time-dependent LLG equation of motion with inertia term.
%
% This file contains the logic and main calculations. Here, the full Green's function is used.

%Integration values in time and energy
t=[0:tstep:tmax];
step=0.1;
w=[-20:step:20];

%Defining a vector for the spin
Sx=zeros(1,length(t));
Sy=zeros(1,length(t));
Sz=zeros(1,length(t));
Sx(1)=-Sxy*sin(wL*t0);
Sy(1)=Sxy*cos(wL*t0);
Sz(1)=Sz0;

%Derivative used in tdLLGinertia calculations
dSx=-Sxy*wL*cos(wL*t0);
dSy=-Sxy*wL*sin(wL*t0);
dSz=0;
ddSx=Sxy*wL^2*sin(wL*t0);
ddSy=-Sxy*wL^2*cos(wL*t0);
ddSz=0;

%Polarized gamma for leads with different magnetization
polarizedgamma

%Defining fermi-functions
fermifunction

%Initiating some variables for speed
initiatingvariables

%Run ODE-solver with Heuns method, i.e., S(n+1) = S(n) + h/2*(dS(n,S(n))+dS(n+h,S(n)+dS(n,S(n)))).
%First loop over each timestep in t.
for j=1:length(t)

    %For calculations before t0 is performed with no exchange coupling J
    if t(j)>=t0
        J=J0;
    else
        J=0;
    end

    %Time the calculation for each timestep
    if(t(j)==floor(t(j)))
        timestep=t(j);
        timeused=toc;
        disp(['Timestep: ' num2str(timestep)])
        disp(['Timeused: ' num2str(timeused)])
    end

    %Tau indicates integration from minus infinity to t. Here with a
    %cut-off at tback.
    tau = [-tstep2*tback+t(j):tstep2:t(j)];

    %Initiate first spin in Heuns method
    S = [Sx(j), Sy(j), Sz(j)];
    dS = [dSx, dSy, dSz];
    ddS = [ddSx, ddSy, ddSz];

    for i=1:tback+1
        %Green's function of (t,tau) for each timestep tau with integration over energies w/omega
        [G0less(i), G0great(i), G1xless(i), G1xgreat(i), G1yless(i), G1ygreat(i), G1zless(i), G1zgreat(i)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

        %Green's function of (tau,t) for each timestep tau with integration over energies w/omega
        [G0less2(i), G0great2(i), G1xless2(i), G1xgreat2(i), G1yless2(i), G1ygreat2(i), G1zless2(i), G1zgreat2(i)]= greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

        %Calculate self-energy K
        selfenergyK
    end

    %Calculate charge and spin currents by integration over tau
    currents

    %Calculate the exchange interaction given the Green's functions
    exchangeinteraction

    %Calculate the internal field given the Green's functions
    internalfield

    %Calculate the magnetic occupation of the QD
    [mx, my, mz] = QDmagneticoccupation(t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

    %Calculate the effective fields
    effectivefields

    %Spin equation of motion
    [dSx1, dSy1, dSz1, ddSx1, ddSy1, ddSz1] = spinequationofmotiontdLLGinertia(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, dS, ddS, tau, t(j));

    %Saving the fields to plot
    savingfieldstdLLGintertia

    %Calculate the second spin for Heuns method iteration
    Sx2=Sx(j)+dSx1;
    Sy2=Sy(j)+dSy1;
    Sz2=Sz(j)+dSz1;
    S2 = [Sx2, Sy2, Sz2];
    dS1 = [dSx1, dSy1, dSz1];
    ddS1 = [ddSx1, ddSy1, ddSz1];

    %Second step in Heuns method iteration
    for i=1:tback+1

        %Green's function of (t,tau) for each timestep tau with integration over energies w/omega
        [G0less(i), G0great(i), G1xless(i), G1xgreat(i), G1yless(i), G1ygreat(i), G1zless(i), G1zgreat(i)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, mu, w, fermi, S2, J);

        %Green's function of (tau, t) for each timestep tau with integration over energies w/omega
        [G0less2(i), G0great2(i), G1xless2(i), G1xgreat2(i), G1yless2(i), G1ygreat2(i), G1zless2(i), G1zgreat2(i)]= greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S2, J);

        %Calculate self-energy K
        selfenergyK

    end

    exchangeinteraction

    internalfield

    [mx, my, mz] = QDmagneticoccupation(t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S2, J);

    effectivefields

    %Spin equation of motion
    [dSx2, dSy2, dSz2, ddSx2, ddSy2, ddSz2] = spinequationofmotiontdLLGinertia(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S2, dS1, ddS1, tau, t(j));

    dSx=(dSx1+dSx2)/2;
    dSy=(dSy1+dSy2)/2;
    dSz=(dSz1+dSz2)/2;

    ddSx=(ddSx1+ddSx2)/2;
    ddSy=(ddSy1+ddSy2)/2;
    ddSz=(ddSz1+ddSz2)/2;

    %Calculate final spin and normalize it
    normalizingspin
end

%Converting time and currents units
currentSIconvert
