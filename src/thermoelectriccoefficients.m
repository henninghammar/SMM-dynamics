function [T_up, T_down, GL, GR, GSL, GSR, ScL, ScR, SsL, SsR, kL, kR, ZTL, ZTR, ZTSL, ZTSR] = thermoelectriccoefficients(pL, pR, gamma, mu, eps, w, dw, J, S, wL, beta, T)
    [G0ret, G0adv, G1xret, G1xadv, G1yret, G1yadv, G1zret, G1zadv] = stationaryretardedgreensfunction(pL, pR, gamma, mu, eps, w, J, S, wL, beta);

    g0L=gamma; %Gamma left
    g0R=gamma; %Gamma righ
    gSL=pL*g0L;
    gSR=pR*g0R;
    gR_up=g0R/2*(1+pR); %Gamma up
    gR_down=g0R/2*(1-pR); %Gamma down
    gL_up=g0L/2*(1+pL); %Gamma up
    gL_down=g0L/2*(1-pL); %Gamma down

    dfermiL = -exp(beta(1).*(w+mu(1)))./(1+exp(beta(1).*(w+mu(1)))).^2;
    dfermiR = -exp(beta(2).*(w+mu(2)))./(1+exp(beta(2).*(w+mu(2)))).^2;

    Gret_up = G0ret + G1zret;
    Gret_down = G0ret - G1zret;
    Gadv_up = G0adv + G1zadv;
    Gadv_down = G0adv - G1zadv;

    %Transmission
    T_up = Gadv_up.*gR_up.*Gret_up.*gL_up;
    T_down = Gadv_down.*gR_down.*Gret_down.*gL_down;

    %Kinetic coefficients
    L0_upL = -trapz(w, 1/(2*pi).*dfermiL.*T_up);
    L1_upL = -trapz(w, 1/(2*pi).*dfermiL.*(w-mu(1)).*T_up);
    L2_upL = -trapz(w, 1/(2*pi).*dfermiL.*(w-mu(1)).^2.*T_up);
    L0_upR = -trapz(w, 1/(2*pi).*dfermiR.*T_up);
    L1_upR = -trapz(w, 1/(2*pi).*dfermiR.*(w-mu(2)).*T_up);
    L2_upR = -trapz(w, 1/(2*pi).*dfermiR.*(w-mu(2)).^2.*T_up);
    L0_downL = -trapz(w, 1/(2*pi).*dfermiL.*T_down);
    L1_downL = -trapz(w, 1/(2*pi).*dfermiL.*(w-mu(1)).*T_down);
    L2_downL = -trapz(w, 1/(2*pi).*dfermiL.*(w-mu(1)).^2.*T_down);
    L0_downR = -trapz(w, 1/(2*pi).*dfermiR.*T_down);
    L1_downR = -trapz(w, 1/(2*pi).*dfermiR.*(w-mu(2)).*T_down);
    L2_downR = -trapz(w, 1/(2*pi).*dfermiR.*(w-mu(2)).^2.*T_down);

    %Conductance
    GL = L0_upL + L0_downL;
    GR = L0_upR + L0_downR;
    GSL = (L0_upL - L0_downL)/2;
    GSR = (L0_upR - L0_downR)/2;

    %Seebeck coefficients
    S_upL = - 1/T(1).*L1_upL/L0_upL;
    S_downL = - 1/T(1).*L1_downL/L0_downL;
    S_upR = - 1/T(2).*L1_upR/L0_upR;
    S_downR = - 1/T(2).*L1_downR/L0_downR;

    %Charge and spin Seebeck
    ScL = (S_upL + S_downL)/2;
    ScR = (S_upR + S_downR)/2;
    SsL = (S_upL - S_downL)/2;
    SsR = (S_upR - S_downR)/2;

    %Electronic thermal Conductance
    kL = 1/T(1)*(L2_upL-L1_upL^2/L0_upL+L2_downL-L1_downL^2/L0_downL);
    kR = 1/T(1)*(L2_upR-L1_upR^2/L0_upR+L2_downR-L1_downR^2/L0_downR);

    %Figure of merit
    ZTL = abs(GL)*ScL^2*T(1)/kL;
    ZTR = abs(GR)*ScR^2*T(2)/kR;
    ZTSL = abs(GSL)*SsL^2*T(1)/kL;
    ZTSR = abs(GSR)*SsR^2*T(2)/kR;
end
