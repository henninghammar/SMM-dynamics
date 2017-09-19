% Non-polarized single-molecule magnet solved with constant LLG parameters
% ================================================
% Author: Henning Hammar
% ------------------
% Calculation of the spin dynamics and current through a degenerate single-molecule magnet
% solved with constant LLG parameters.
%
% This file contains the logic and main calculations. Here, the non-polarized Green's function is used.

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

%Derivative used in LLG calculations
dSx=-Sxy*wL*cos(wL*t0);
dSy=-Sxy*wL*sin(wL*t0);
dSz=0;

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

    S = [Sx(j), Sy(j), Sz(j)];
    dS = [dSx, dSy, dSz];

    %Calculate the interactionparameters for the first timestep
    if j == 1
        [jH, jIxx, jIyy, jIzz, jIxy, jIxz, jIyz, jDMx, jDMy, jDMz, ejx, ejy, ejz, GjH, GjIxx, GjIyy, GjIzz, GjIxy, GjIxz, GjIyz, GjDMx, GjDMy, GjDMz] = interactionparameters(pL, pR, gamma, eV, eps, epsilon, w, step, J, S, wL, beta);

        %Calculate the magnetic occupation of the QD
        [mx, my, mz] = degenerateQDmagneticoccupation(w, epsilon, g0tot, J, fermi, S);

        %Calculate the effective fields
        effectivefields

        jIyx = jIxy;
        jIzx = jIxz;
        jIzy = jIyz;
        GjIyx = GjIxy;
        GjIzx = GjIxz;
        GjIzy = GjIyz;
    end

    %Spin equation of motion
    [dSx1, dSy1, dSz1] = spinequationofmotionLLGconstant(Beffx, Beffy, Beffz, jH, jIxx, jIyy, jIzz, jIxy, jIxz, jIyz, jIyx, jIzx, jIzy, jDMx, jDMy, jDMz, GjH, GjIxx, GjIyy, GjIzz, GjIxy, GjIxz, GjIyx, GjIzx, GjIzy, GjIyz, GjDMx, GjDMy, GjDMz, S, dS);

    %Calculate the second spin for Heuns method iteration
    Sx2=Sx(j)+dSx1;
    Sy2=Sy(j)+dSy1;
    Sz2=Sz(j)+dSz1;
    S2 = [Sx2, Sy2, Sz2];
    dS1 = [dSx1, dSy1, dSz1];

    %Spin equation of motion
    [dSx2, dSy2, dSz2] = spinequationofmotionLLGconstant(Beffx, Beffy, Beffz, jH, jIxx, jIyy, jIzz, jIxy, jIxz, jIyz, jIyx, jIzx, jIzy, jDMx, jDMy, jDMz, GjH, GjIxx, GjIyy, GjIzz, GjIxy, GjIxz, GjIyx, GjIzx, GjIzy, GjIyz, GjDMx, GjDMy, GjDMz, S2, dS1);

    dSx=(dSx1+dSx2)/2;
    dSy=(dSy1+dSy2)/2;
    dSz=(dSz1+dSz2)/2;

    %Calculate final spin and normalize it
    normalizingspin
end

%Converting time and currents to SI-units
SIconvert
