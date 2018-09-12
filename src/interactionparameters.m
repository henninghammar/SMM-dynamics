function [JH, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Dx, Dy, Dz, ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, mu, eps, epsilon, w, dw, J, S, wL, beta)
    [G0less, G0great, G1xless, G1xgreat, G1yless, G1ygreat, G1zless, G1zgreat] = stationarygreensfunction(pL, pR, gamma, mu, eps, w, J, S, wL, beta);

    delta = 10^-5;
    %
    % %Calculation method 1
    % for i = 1:length(w)
    %     w2 = w(i);
    %     W(i,:) = 1./(w-w2+1i*delta);
    %     Wz(i,:) = 1./(w-w2+delta);
    %     WG(i,:) = 1./(w-w2+1i*delta).^2;
    % end
    %
    % ejx=-J/(2*pi^2).*(...
    %     (epsilon*G0less+wL/2.*G1zless)*W*G1xgreat'-(epsilon*G0great+wL/2.*G1zgreat)*W*G1xless'...
    %     +(epsilon*G1xless-1i.*G1yless)*W*G0great'-(epsilon*G1xgreat-1i.*G1ygreat)*W*G0less'...
    %     -1i.*(epsilon*G1yless+1i.*G1xless)*W*G1zgreat'+1i.*(epsilon*G1zless+wL/2)*W*G1ygreat'...
    %     +1i.*(epsilon*G1ygreat+1i.*G1xgreat)*W*G1zless'-1i.*(epsilon*G1zgreat+wL/2)*W*G1yless');
    % ejy=-J/(2*pi^2).*(...
    %     (epsilon*G0less+wL/2.*G1zless)*W*G1ygreat'-(epsilon*G0great+wL/2.*G1zgreat)*W*G1yless'...
    %     +(epsilon*G1yless+1i.*G1xless)*W*G0great'-(epsilon*G1ygreat+1i.*G1xgreat)*W*G0less'...
    %     -1i.*(epsilon*G1zless+wL/2)*W*G1xgreat'+1i.*(epsilon*G1xless-1i.*G1yless)*W*G1zgreat'...
    %     +1i.*(epsilon*G1zgreat+wL/2)*W*G1xless'-1i.*(epsilon*G1xgreat-1i.*G1ygreat)*W*G1zless');
    % ejz=-J/(2*pi^2).*(...
    %     (epsilon*G0less+wL/2.*G1zless)*Wz*G1zgreat'-(epsilon*G0great+wL/2.*G1zgreat)*Wz*G1zless'...
    %     +(epsilon*G1zless+wL/2)*Wz*G0great'-(epsilon*G1zgreat+wL/2)*Wz*G0less'...
    %     -1i.*(epsilon*G1xless-1i.*G1yless)*Wz*G1ygreat'+1i.*(epsilon*G1yless+1i.*G1xless)*Wz*G1xgreat'...
    %     +1i.*(epsilon*G1xgreat-1i.*G1ygreat)*Wz*G1yless'-1i.*(epsilon*G1ygreat+1i.*G1xgreat)*Wz*G1xless');
    %
    % ejx = real(ejx)*dw^2;
    % ejy = real(ejy)*dw^2;
    % ejz = real(ejz)*dw^2;
    %
    % JH = -J^2/(4*pi^2).*(G0less*W*G0great'-G0great*W*G0less'-G1xless*W*G1xgreat'+G1xgreat*W*G1xless'...
    %   -G1yless*W*G1ygreat'+G1ygreat*W*G1yless'-G1zless*W*G1zgreat'+G1zgreat*W*G1zless');
    % Ixx = -J^2/(2*pi^2).*(G1xless*W*G1xgreat'-G1xgreat*W*G1xless');
    % Iyy = -J^2/(2*pi^2).*(G1yless*W*G1ygreat'-G1ygreat*W*G1yless');
    % Izz = -J^2/(2*pi^2).*(G1zless*W*G1zgreat'-G1zgreat*W*G1zless');
    % Ixz = -J^2/(4*pi^2).*(G1xless*W*G1zgreat'-G1xgreat*W*G1zless'+G1zless*W*G1xgreat'-G1zgreat*W*G1xless');
    % Ixy = -J^2/(4*pi^2).*(G1xless*W*G1ygreat'-G1xgreat*W*G1yless'+G1yless*W*G1xgreat'-G1ygreat*W*G1xless');
    % Iyz = -J^2/(4*pi^2).*(G1yless*W*G1zgreat'-G1ygreat*W*G1zless'+G1zless*W*G1ygreat'-G1zgreat*W*G1yless');
    % Dx = 1i*J^2/(4*pi^2).*(G0less*W*G1xgreat'-G0great*W*G1xless'-G1xless*W*G0great'+G1xgreat*W*G0less');
    % Dy = 1i*J^2/(4*pi^2).*(G0less*W*G1ygreat'-G0great*W*G1yless'-G1yless*W*G0great'+G1ygreat*W*G0less');
    % Dz = 1i*J^2/(4*pi^2).*(G0less*W*G1zgreat'-G0great*W*G1zless'-G1zless*W*G0great'+G1zgreat*W*G0less');
    %
    % JH = real(JH)*dw^2;
    % Ixx = real(Ixx)*dw^2;
    % Iyy = real(Iyy)*dw^2;
    % Izz = real(Izz)*dw^2;
    % Ixz = real(Ixz)*dw^2;
    % Ixy = real(Ixy)*dw^2;
    % Iyz = real(Iyz)*dw^2;
    % Dx = real(Dx)*dw^2;
    % Dy = real(Dy)*dw^2;
    % Dz = real(Dz)*dw^2;
    %
    % GJH = -J^2/(4*pi^2).*(G0less*WG*G0great'-G0great*WG*G0less'-G1xless*WG*G1xgreat'+G1xgreat*WG*G1xless'...
    %   -G1yless*WG*G1ygreat'+G1ygreat*WG*G1yless'-G1zless*WG*G1zgreat'+G1zgreat*WG*G1zless');
    % GIxx = -J^2/(2*pi^2).*(G1xless*WG*G1xgreat'-G1xgreat*WG*G1xless');
    % GIyy = -J^2/(2*pi^2).*(G1yless*WG*G1ygreat'-G1ygreat*WG*G1yless');
    % GIzz = -J^2/(2*pi^2).*(G1zless*WG*G1zgreat'-G1zgreat*WG*G1zless');
    % GIxz = -J^2/(4*pi^2).*(G1xless*WG*G1zgreat'-G1xgreat*WG*G1zless'+G1zless*WG*G1xgreat'-G1zgreat*WG*G1xless');
    % GIxy = -J^2/(4*pi^2).*(G1xless*WG*G1ygreat'-G1xgreat*WG*G1yless'+G1yless*WG*G1xgreat'-G1ygreat*WG*G1xless');
    % GIyz = -J^2/(4*pi^2).*(G1yless*WG*G1zgreat'-G1ygreat*WG*G1zless'+G1zless*WG*G1ygreat'-G1zgreat*WG*G1yless');
    % GDx = 1i*J^2/(4*pi^2).*(G0less*WG*G1xgreat'-G0great*WG*G1xless'-G1xless*WG*G0great'+G1xgreat*WG*G0less');
    % GDy = 1i*J^2/(4*pi^2).*(G0less*WG*G1ygreat'-G0great*WG*G1yless'-G1yless*WG*G0great'+G1ygreat*WG*G0less');
    % GDz = 1i*J^2/(4*pi^2).*(G0less*WG*G1zgreat'-G0great*WG*G1zless'-G1zless*WG*G0great'+G1zgreat*WG*G0less');
    %
    % GJH = -imag(GJH)*dw^2;
    % GIxx = -imag(GIxx)*dw^2;
    % GIyy = -imag(GIyy)*dw^2;
    % GIzz = -imag(GIzz)*dw^2;
    % GIxz = -imag(GIxz)*dw^2;
    % GIxy = -imag(GIxy)*dw^2;
    % GIyz = -imag(GIyz)*dw^2;
    % GDx = -imag(GDx)*dw^2;
    % GDy = -imag(GDy)*dw^2;
    % GDz = -imag(GDz)*dw^2;

    %Calculation method 2
    w2 = linspace(-10,10,203);
    for i = 1:length(w2)
        [G0less2, G0great2, G1xless2, G1xgreat2, G1yless2, G1ygreat2, G1zless2, G1zgreat2] = stationarygreensfunction(pL, pR, gamma, mu, eps, w2(i), J, S, wL, beta);

        ejxbare=-J^2/(2*pi^2).*((epsilon.*G0less+wL/2.*G1zless).*G1xgreat2-(epsilon.*G0great+wL/2.*G1zgreat).*G1xless2...
            +(epsilon.*G1xless-1i.*G1yless).*G0great2-(epsilon.*G1xgreat-1i.*G1ygreat).*G0less2...
            -1i.*(epsilon.*G1yless+1i.*G1xless).*G1zgreat2+1i.*(epsilon.*G1zless+wL/2).*G1ygreat2...
            +1i.*(epsilon.*G1ygreat+1i.*G1xgreat).*G1zless2-1i.*(epsilon.*G1zgreat+wL/2).*G1yless2)./(w-w2(i)+1i*delta);
        ejybare=-J^2/(2*pi^2).*((epsilon*G0less+wL/2.*G1zless).*G1ygreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1yless2...
            +(epsilon*G1yless+1i.*G1xless).*G0great2-(epsilon*G1ygreat+1i.*G1xgreat).*G0less2...
            -1i.*(epsilon*G1zless+wL/2).*G1xgreat2+1i.*(epsilon*G1xless-1i.*G1yless).*G1zgreat2...
            +1i.*(epsilon*G1zgreat+wL/2).*G1xless2-1i.*(epsilon*G1xgreat-1i.*G1ygreat).*G1zless2)./(w-w2(i)+1i*delta);
        ejzbare=-J^2/(2*pi^2).*((epsilon*G0less+wL/2.*G1zless).*G1zgreat2-(epsilon*G0great+wL/2.*G1zgreat).*G1zless2...
            +(epsilon*G1zless+wL/2).*G0great2-(epsilon*G1zgreat+wL/2).*G0less2...
            -1i.*(epsilon*G1xless-1i.*G1yless).*G1ygreat2+1i.*(epsilon*G1yless+1i.*G1xless).*G1xgreat2...
            +1i.*(epsilon*G1xgreat-1i.*G1ygreat).*G1yless2-1i.*(epsilon*G1ygreat+1i.*G1xgreat).*G1xless2)./(w-w2(i)+delta);

        ejxvect(i) = trapz(w, ejxbare);
        ejyvect(i) = trapz(w, ejybare);
        ejzvect(i) = trapz(w, ejzbare);

        JHbare = -J^2/(4*pi^2).*(G0less.*G0great2-G0great.*G0less2-G1xless.*G1xgreat2+G1xgreat.*G1xless2...
          -G1yless.*G1ygreat2+G1ygreat.*G1yless2-G1zless.*G1zgreat2+G1zgreat.*G1zless2)./(w-w2(i)+1i*delta);
        Ixxbare = -J^2/(2*pi^2).*(G1xless.*G1xgreat2-G1xgreat.*G1xless2)./(w-w2(i)+1i*delta);
        Iyybare = -J^2/(2*pi^2).*(G1yless.*G1ygreat2-G1ygreat.*G1yless2)./(w-w2(i)+1i*delta);
        Izzbare = -J^2/(2*pi^2).*(G1zless.*G1zgreat2-G1zgreat.*G1zless2)./(w-w2(i)+1i*delta);
        Ixzbare = -J^2/(4*pi^2).*(G1xless.*G1zgreat2-G1xgreat.*G1zless2+G1zless.*G1xgreat2-G1zgreat.*G1xless2)./(w-w2(i)+1i*delta);
        Ixybare = -J^2/(4*pi^2).*(G1xless.*G1ygreat2-G1xgreat.*G1yless2+G1yless.*G1xgreat2-G1ygreat.*G1xless2)./(w-w2(i)+1i*delta);
        Iyzbare = -J^2/(4*pi^2).*(G1yless.*G1zgreat2-G1ygreat.*G1zless2+G1zless.*G1ygreat2-G1zgreat.*G1yless2)./(w-w2(i)+1i*delta);
        Dxbare = -1i*J^2/(4*pi^2).*(G0less.*G1xgreat2-G0great.*G1xless2-G1xless.*G0great2+G1xgreat.*G0less2)./(w-w2(i)+1i*delta);
        Dybare = -1i*J^2/(4*pi^2).*(G0less.*G1ygreat2-G0great.*G1yless2-G1yless.*G0great2+G1ygreat.*G0less2)./(w-w2(i)+1i*delta);
        Dzbare = -1i*J^2/(4*pi^2).*(G0less.*G1zgreat2-G0great.*G1zless2-G1zless.*G0great2+G1zgreat.*G0less2)./(w-w2(i)+1i*delta);

        JHvect(i) = trapz(w, JHbare);
        Ixxvect(i) = trapz(w, Ixxbare);
        Iyyvect(i) = trapz(w, Iyybare);
        Izzvect(i) = trapz(w, Izzbare);
        Ixzvect(i) = trapz(w, Ixzbare);
        Ixyvect(i) = trapz(w, Ixybare);
        Iyzvect(i) = trapz(w, Iyzbare);
        Dxvect(i) = trapz(w, Dxbare);
        Dyvect(i) = trapz(w, Dybare);
        Dzvect(i) = trapz(w, Dzbare);

        GJHbare = -J^2/(4*pi^2).*(G0less.*G0great2-G0great.*G0less2-G1xless.*G1xgreat2+G1xgreat.*G1xless2...
          -G1yless.*G1ygreat2+G1ygreat.*G1yless2-G1zless.*G1zgreat2+G1zgreat.*G1zless2)./(w-w2(i)+1i*delta).^2;
        GIxxbare = -J^2/(2*pi^2).*(G1xless.*G1xgreat2-G1xgreat.*G1xless2)./(w-w2(i)+1i*delta).^2;
        GIyybare = -J^2/(2*pi^2).*(G1yless.*G1ygreat2-G1ygreat.*G1yless2)./(w-w2(i)+1i*delta).^2;
        GIzzbare = -J^2/(2*pi^2).*(G1zless.*G1zgreat2-G1zgreat.*G1zless2)./(w-w2(i)+1i*delta).^2;
        GIxzbare = -J^2/(4*pi^2).*(G1xless.*G1zgreat2-G1xgreat.*G1zless2+G1zless.*G1xgreat2-G1zgreat.*G1xless2)./(w-w2(i)+1i*delta).^2;
        GIxybare = -J^2/(4*pi^2).*(G1xless.*G1ygreat2-G1xgreat.*G1yless2+G1yless.*G1xgreat2-G1ygreat.*G1xless2)./(w-w2(i)+1i*delta).^2;
        GIyzbare = -J^2/(4*pi^2).*(G1yless.*G1zgreat2-G1ygreat.*G1zless2+G1zless.*G1ygreat2-G1zgreat.*G1yless2)./(w-w2(i)+1i*delta).^2;
        GDxbare = -1i*J^2/(4*pi^2).*(G0less.*G1xgreat2-G0great.*G1xless2-G1xless.*G0great2+G1xgreat.*G0less2)./(w-w2(i)+1i*delta).^2;
        GDybare = -1i*J^2/(4*pi^2).*(G0less.*G1ygreat2-G0great.*G1yless2-G1yless.*G0great2+G1ygreat.*G0less2)./(w-w2(i)+1i*delta).^2;
        GDzbare = -1i*J^2/(4*pi^2).*(G0less.*G1zgreat2-G0great.*G1zless2-G1zless.*G0great2+G1zgreat.*G0less2)./(w-w2(i)+1i*delta).^2;

        GJHvect(i) = -imag(trapz(w, GJHbare));
        GIxxvect(i) = -imag(trapz(w, GIxxbare));
        GIyyvect(i) = -imag(trapz(w, GIyybare));
        GIzzvect(i) = -imag(trapz(w, GIzzbare));
        GIxzvect(i) = -imag(trapz(w, GIxzbare));
        GIxyvect(i) = -imag(trapz(w, GIxybare));
        GIyzvect(i) = -imag(trapz(w, GIyzbare));
        GDxvect(i) = -imag(trapz(w, GDxbare));
        GDyvect(i) = -imag(trapz(w, GDybare));
        GDzvect(i) = -imag(trapz(w, GDzbare));
    end

    ejx = trapz(w2,ejxvect);
    ejy = trapz(w2,ejyvect);
    ejz = trapz(w2,ejzvect);
    JH = real(trapz(w2, JHvect));
    Ixx = real(trapz(w2, Ixxvect));
    Iyy = real(trapz(w2, Iyyvect));
    Izz = real(trapz(w2, Izzvect));
    Ixz = real(trapz(w2, Ixzvect));
    Ixy = real(trapz(w2, Ixyvect));
    Iyz = real(trapz(w2, Iyzvect));
    Dx = real(trapz(w2, Dxvect));
    Dy = real(trapz(w2, Dyvect));
    Dz = real(trapz(w2, Dzvect));
    GJH = trapz(w2, GJHvect);
    GIxx = trapz(w2, GIxxvect);
    GIyy = trapz(w2, GIyyvect);
    GIzz = trapz(w2, GIzzvect);
    GIxz = trapz(w2, GIxzvect);
    GIxy = trapz(w2, GIxyvect);
    GIyz = trapz(w2, GIyzvect);
    GDx = trapz(w2, GDxvect);
    GDy = trapz(w2, GDyvect);
    GDz = trapz(w2, GDzvect);
end
