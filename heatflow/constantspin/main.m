% Pulse with a constant SMM
% ================================================
% Author: Henning Hammar
% ------------------
% Calculation of the current and other observables for a pulse with a constant SMM.
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
S = [Sx(1), Sy(1), Sz(1)];

%Polarized gamma for leads with different magnetization
polarizedgamma

%Defining fermi-functions
fermifunction

%Initiating some variables for speed
initiatingvariables

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

    %Calculate the occupations and DOS of the QD
    [n(j), mx, my, mz, DOS(j)] = QDoccupation(t(j), t0, t1, eps, g, g0, gS, mu, w, fermi, S, J);
    renormalization

    %Calculate charge and spin currents by integration over tau
    currents
    energycurrents

    %Calculate the exchange interaction given the Green's functions
    exchangeinteraction

    %Calculate the internal field given the Green's functions
    internalfield

    %Calculate the effective fields
    effectivefields

    %Saving the fields to plot
    %savingfields
    %The fields acting on the spin
    ejvect(:,j)=[ejx, ejy, ejz];

    mvect(:,j)=[mx, my, mz];

    Beffvect(:,j)=[Beffx, Beffy, Beffz];

    jHt(j)=trapz(tau,jH);
    jDMvect(:,j)=[trapz(tau,jDMx), trapz(tau,jDMy), trapz(tau,jDMz)];
    jIvect(:,j)=[trapz(tau,jIxx), trapz(tau,jIyy), trapz(tau,jIzz), trapz(tau,jIxy), trapz(tau,jIxz), trapz(tau,jIyz)];

end

%Converting currents units
currentSIconvert
