%The core program

%Integration values in time and energy
tstep2=0.1;
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
    
    for i=1:tback+1 
        
%         S = [Sx(j), Sy(j), Sz(j)];
%         
%         %Green's function of (t',t)
%         [G0less(j), G0great(j), G1xless(j), G1xgreat(j), G1yless(j), G1ygreat(j), G1zless(j), G1zgreat(j)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, eV, w, fermi, S, J);
%         
%         %Green's function of (t',t)
%         [G0less2(j), G0great2(j), G1xless2(j), G1xgreat2(j), G1yless2(j), G1ygreat2(j), G1zless2(j), G1zgreat2(j)] = greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, eV, w, fermi, S, J);
        
        for k=1:2
            
            %Green's function of (t',t)
            baregreensfunction1
            fullgreensfunction1

            %Green's function of (t,t')
            baregreensfunction2
            fullgreensfunction2
        end
        
        greensfunctionenergyintegration
        
        %Calculate self-energy K
        selfenergyK
    end
    
    %Calculate charge and spin currents
    currents
    
    %Calculate the exchange interaction given the Green's functions
    exchangeinteraction
    
    %Calculate the internal field given the Green's functions
    internalfield
    
    %Spin equation of motion
    spinequationofmotion
    
    %Saving the fields to plot
    savingfields
    
    %Calculate new spin for Heuns method iteration
    Sx2=Sx(j)+dSx1;
    Sy2=Sy(j)+dSy1;
    Sz2=Sz(j)+dSz1;
    
    %Second step in Heuns method iteration
    for i=1:tback+1 

        %S2 = [Sx2, Sy2, Sz2];
        
        %Green's function of (t',t)
        %[G0less(j), G0great(j), G1xless(j), G1xgreat(j), G1yless(j), G1ygreat(j), G1zless(j), G1zgreat(j)] = greensfunction(t(j), tau(i), t0, t1, eps, g, g0, gS, eV, w, fermi, S2, J);
        
        %Green's function of (t',t)
        %[G0less2(j), G0great2(j), G1xless2(j), G1xgreat2(j), G1yless2(j), G1ygreat2(j), G1zless2(j), G1zgreat2(j)] = greensfunction(tau(i), t(j), t0, t1, eps, g, g0, gS, eV, w, fermi, S2, J);
        
        for k=1:2
            
            %Green's function of (t',t)
            baregreensfunction1
            fullgreensfunction1alt2
            
            selfenergyK
            
            %Green's function of (t,t')
            baregreensfunction2
            fullgreensfunction2alt2
        end
        
        greensfunctionenergyintegration
        
        %Calculate self-energy K
        selfenergyK
    
    end
    
    exchangeinteraction
    
    internalfield
    
    spinequationofmotionalt2
    
    %Calculate final spin and normalize it
    normalizingspin
end

%Converting time and currents to SI-units
SIconvert