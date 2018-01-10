% Polarized single-molecule magnet with pulse
% ================================================
% Author: Henning Hammar
% ------------------
% Old setup before rewriting the document.
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

%Polarized gamma for leads with different magnetization
polarizedgamma

%Defining fermi-functions
fermifunction

%Initiating some variables for speed
initiatingvariables

%Double gamma as in old code
g = g*2;

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

    for i=1:tback+1
        %Green's function of (t,tau) for each timestep tau with integration over energies w/omega
        [G0less(i), G0great(i), G1xless(i), G1xgreat(i), G1yless(i), G1ygreat(i), G1zless(i), G1zgreat(i)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, eV, w, fermi, S, J);

        %Green's function of (tau,t) for each timestep tau with integration over energies w/omega
        [G0less2(i), G0great2(i), G1xless2(i), G1xgreat2(i), G1yless2(i), G1ygreat2(i), G1zless2(i), G1zgreat2(i)] = greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, eV, w, fermi, S, J);

        %Calculate self-energy K
        selfenergyK
    end

    %Calculate charge and spin currents by integration over tau
    currents

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

    ST = [SxT, SyT, SzT];

    %Calculate the exchange interaction given the Green's functions
    exchangeinteraction

    %Calculate the internal field given the Green's functions
    internalfield

    Beffx = ejx;
    Beffy = ejy;
    Beffz = wL + ejz;

    %Spin equation of motion
    [dSx1, dSy1, dSz1] = spinequationofmotion(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S, ST, tau);

    %Saving the fields to plot
    %Saving the interesting fields
    %Fields with spin
    SBx(j)=-Beffz.*Sy(j)+Beffy.*Sz(j);
    SBy(j)=-Beffx.*Sz(j)+Beffz.*Sx(j);
    SBz(j)=-Beffy.*Sx(j)+Beffx.*Sy(j);

    Sejx(j)=ejz.*Sy(j)-ejy.*Sz(j);
    Sejy(j)=ejx.*Sz(j)-ejz.*Sx(j);
    Sejz(j)=ejy.*Sx(j)-ejx.*Sy(j);

    dSxbarejH=jH.*(Sy(j).*SzT-Sz(j).*SyT);
    dSybarejH=jH.*(Sz(j).*SxT-Sx(j).*SzT);
    dSzbarejH=jH.*(Sx(j).*SyT-Sy(j).*SxT);
    SjHx(j)=trapz(tau,dSxbarejH);
    SjHy(j)=trapz(tau,dSybarejH);
    SjHz(j)=trapz(tau,dSzbarejH);

    dSxbarejDMx=-Sy(j).*(jDMx.*SyT-jDMy.*SxT)+Sz(j).*(jDMz.*SxT-jDMx.*SzT);
    dSybarejDMy=-Sz(j).*(jDMy.*SzT-jDMz.*SyT)+Sx(j).*(jDMx.*SyT-jDMy.*SxT);
    dSzbarejDMz=-Sx(j).*(jDMz.*SxT-jDMx.*SzT)+Sy(j).*(jDMy.*SzT-jDMz.*SyT);
    SjDMx(j)=trapz(tau,dSxbarejDMx);
    SjDMy(j)=trapz(tau,dSybarejDMy);
    SjDMz(j)=trapz(tau,dSzbarejDMz);

    dSxbarejIx=Sy(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT)-Sz(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT);
    dSybarejIy=Sz(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT)-Sx(j).*(jIzx.*SxT+jIzy.*SyT+jIzz.*SzT);
    dSzbarejIz=Sx(j).*(jIyx.*SxT+jIyy.*SyT+jIyz.*SzT)-Sy(j).*(jIxx.*SxT+jIxy.*SyT+jIxz.*SzT);
    SjIx(j)=trapz(tau,dSxbarejIx);
    SjIy(j)=trapz(tau,dSybarejIy);
    SjIz(j)=trapz(tau,dSzbarejIz);

    %The fields acting on the spin
    ejvect(:,j)=[ejx, ejy, ejz];

    Beffvect(:,j)=[Beffx, Beffy, Beffz];

    jHt(j)=trapz(tau,jH);
    jDMvect(:,j)=[trapz(tau,jDMx), trapz(tau,jDMy), trapz(tau,jDMz)];
    jIvect(:,j)=[trapz(tau,jIxx), trapz(tau,jIyy), trapz(tau,jIzz), trapz(tau,jIxy), trapz(tau,jIxz), trapz(tau,jIyz)];

    %Calculate the second spin for Heuns method iteration
    Sx2=Sx(j)+dSx1;
    Sy2=Sy(j)+dSy1;
    Sz2=Sz(j)+dSz1;
    SxT2=[SxT(2:end),Sx2];
    SyT2=[SyT(2:end),Sy2];
    SzT2=[SzT(2:end),Sz2];
    S2 = [Sx2, Sy2, Sz2];
    ST2 = [SxT2, SyT2, SzT2];

    %Second step in Heuns method iteration
    for i=1:tback+1
        %Green's function of (t,tau) for each timestep tau with integration over energies w/omega
        [G0less(i), G0great(i), G1xless(i), G1xgreat(i), G1yless(i), G1ygreat(i), G1zless(i), G1zgreat(i)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, eV, w, fermi, S2, J);

        %Green's function of (tau,t) for each timestep tau with integration over energies w/omega
        [G0less2(i), G0great2(i), G1xless2(i), G1xgreat2(i), G1yless2(i), G1ygreat2(i), G1zless2(i), G1zgreat2(i)] = greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, eV, w, fermi, S, J);

        %Calculate self-energy K
        selfenergyK
    end

    exchangeinteraction

    internalfield

    Beffx = ejx;
    Beffy = ejy;
    Beffz = wL + ejz;

    %Spin equation of motion
    [dSx2, dSy2, dSz2] = spinequationofmotion(Beffx, Beffy, Beffz,jH,jDMx,jDMy,jDMz,jIxx,jIyy,jIzz,jIxy,jIyx,jIxz,jIzx,jIyz,jIzy, S2, ST2, tau);

    %Calculate final spin and normalize it
    normalizingspin
end

%Converting currents units
currentSIconvert
