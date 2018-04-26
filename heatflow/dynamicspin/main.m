% Polarized single-molecule magnet with pulse
% ================================================
% Author: Henning Hammar
% ------------------
% Calculation of spin dynamics and charge and energy currents through a polarized single-molecule magnet
% for a pulse.
%
% This file contains the logic and main calculations. Here, the full Green's function is used.

%Integration values in time and energy
t=[0:tstep:tmax];
step=0.1;
w=[-20:step:20];
step2=0.1;
w2=[-75:step2:75];

%Defining a vector for the spin
Sx=zeros(1,length(t));
Sy=zeros(1,length(t));
Sz=zeros(1,length(t));
Sx(1)=-Sxy*sin(wL*t0);
Sy(1)=Sxy*cos(wL*t0);
Sz(1)=Sz0;

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
    timestep=t(j);
    timeused=toc;
    disp(['Timestep: ' num2str(timestep)])
    disp(['Timeused: ' num2str(timeused)])

    %Tau indicates integration from minus infinity to t. Here with a
    %cut-off at tback.
    tau = [-tstep2*tback+t(j):tstep2:t(j)];

    %Initialize spin for integration over the memory kernel
    if j == 1
       SxT=-Sxy*sin(wL*tau);
       SyT=Sxy*cos(wL*tau);
       SzT=Sz0*ones(1,length(tau));
    else
       SxT=[SxT(2:end),Sx(j)];
       SyT=[SyT(2:end),Sy(j)];
       SzT=[SzT(2:end),Sz(j)];
    end

    %Initiate first spin in Heuns method
    S = [Sx(j), Sy(j), Sz(j)];

    for i=1:tback+1
        %Green's function of (t,tau) for each timestep tau with integration over energies w/omega
        [G0less(i), G0great(i), G1xless(i), G1xgreat(i), G1yless(i), G1ygreat(i), G1zless(i), G1zgreat(i)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

        %Green's function of (tau,t) for each timestep tau with integration over energies w/omega
        [G0less2(i), G0great2(i), G1xless2(i), G1xgreat2(i), G1yless2(i), G1ygreat2(i), G1zless2(i), G1zgreat2(i)] = greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

        %Calculate self-energy K
        selfenergyK
    end

    tauenergy = [-tstepenergy*tbackenergy+t(j):tstepenergy:t(j)];
    for i=1:length(tauenergy)
      [energyKLless(i), energyKLgreat(i), energyKRless(i), energyKRgreat(i)] = energyselfenergyK(w2, t(j), tauenergy(i), t0, t1, mu, T, kB);
      [energyG0less(i), energyG0great(i), energyG1xless(i), energyG1xgreat(i), energyG1yless(i), energyG1ygreat(i), energyG1zless(i), energyG1zgreat(i)] = greensfunction(t(j), tauenergy(i), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);
    end

    %Calculate charge and spin currents by integration over tau
    currents
    energycurrents

    %Calculate the exchange interaction given the Green's functions
    exchangeinteraction

    %Calculate the internal field given the Green's functions
    internalfield

    %Calculate the magnetic occupation of the QD
    [mx, my, mz] = QDmagneticoccupation(t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

    n(j) = QDoccupation(t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);

    %Calculate the effective fields
    effectivefields

    %Spin equation of motion
    [dSx1, dSy1, dSz1] = spinequationofmotion(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, SxT, SyT, SzT, tau);

    %Saving the fields to plot
    savingfields

    %Calculate the second spin for Heuns method iteration
    Sx2=Sx(j)+dSx1;
    Sy2=Sy(j)+dSy1;
    Sz2=Sz(j)+dSz1;
    SxT2=[SxT(2:end),Sx2];
    SyT2=[SyT(2:end),Sy2];
    SzT2=[SzT(2:end),Sz2];
    S2 = [Sx2, Sy2, Sz2];

    %Second step in Heuns method iteration
    for i=1:tback+1

        %Green's function of (t,tau) for each timestep tau with integration over energies w/omega
        [G0less(i), G0great(i), G1xless(i), G1xgreat(i), G1yless(i), G1ygreat(i), G1zless(i), G1zgreat(i)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, mu, w, fermi, S2, J);

        %Green's function of (tau,t) for each timestep tau with integration over energies w/omega
        [G0less2(i), G0great2(i), G1xless2(i), G1xgreat2(i), G1yless2(i), G1ygreat2(i), G1zless2(i), G1zgreat2(i)] = greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S2, J);
    end

    exchangeinteraction

    internalfield

    [mx, my, mz] = QDmagneticoccupation(t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S2, J);

    effectivefields

    %Spin equation of motion
    [dSx2, dSy2, dSz2] = spinequationofmotion(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, SxT2, SyT2, SzT2, tau);

    %Calculate final spin and normalize it
    normalizingspin
end

%Converting currents units
currentSIconvert
