function [JH, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Dx, Dy, Dz, ejx, ejy, ejz, GJH, GIxx, GIyy, GIzz, GIxy, GIxz, GIyz, GDx, GDy, GDz] = interactionparameters(pL, pR, gamma, eV, eps, epsilon, w, dw, J, S, wL, beta)
    stationarygreensfunction

    delta = 10^-5;
    for i = 1:length(w)
        w2 = w(i);
        W(i,:) = 1./(w-w2+1i*delta);
        Wz(i,:) = 1./(w-w2+delta);
        WG(i,:) = 1./(w-w2+1i*delta).^2;
    end

    %W2 = 1./(w-w2+1i*delta).^2;

    ejx=-J/(2*pi^2).*(...
        (epsilon*G0less+wL/2.*G1zless)*W*G1xgreat'-(epsilon*G0great+wL/2.*G1zgreat)*W*G1xless'...
        +(epsilon*G1xless-1i.*G1yless)*W*G0great'-(epsilon*G1xgreat-1i.*G1ygreat)*W*G0less'...
        -1i.*(epsilon*G1yless+1i.*G1xless)*W*G1zgreat'+1i.*(epsilon*G1zless+wL/2)*W*G1ygreat'...
        +1i.*(epsilon*G1ygreat+1i.*G1xgreat)*W*G1zless'-1i.*(epsilon*G1zgreat+wL/2)*W*G1yless');
    ejx = real(ejx)*dw^2;
    ejy=-J/(2*pi^2).*(...
        (epsilon*G0less+wL/2.*G1zless)*W*G1ygreat'-(epsilon*G0great+wL/2.*G1zgreat)*W*G1yless'...
        +(epsilon*G1yless+1i.*G1xless)*W*G0great'-(epsilon*G1ygreat+1i.*G1xgreat)*W*G0less'...
        -1i.*(epsilon*G1zless+wL/2)*W*G1xgreat'+1i.*(epsilon*G1xless-1i.*G1yless)*W*G1zgreat'...
        +1i.*(epsilon*G1zgreat+wL/2)*W*G1xless'-1i.*(epsilon*G1xgreat-1i.*G1ygreat)*W*G1zless');
    ejy = real(ejy)*dw^2;
    ejz=-J/(2*pi^2).*(...
        (epsilon*G0less+wL/2.*G1zless)*Wz*G1zgreat'-(epsilon*G0great+wL/2.*G1zgreat)*Wz*G1zless'...
        +(epsilon*G1zless+wL/2)*Wz*G0great'-(epsilon*G1zgreat+wL/2)*Wz*G0less'...
        -1i.*(epsilon*G1xless-1i.*G1yless)*Wz*G1ygreat'+1i.*(epsilon*G1yless+1i.*G1xless)*Wz*G1xgreat'...
        +1i.*(epsilon*G1xgreat-1i.*G1ygreat)*Wz*G1yless'-1i.*(epsilon*G1ygreat+1i.*G1xgreat)*Wz*G1xless');
    ejz = real(ejz)*dw^2;

    JH = -J^2/(4*pi^2).*(G0less*W*G0great'-G0great*W*G0less'-G1xless*W*G1xgreat'+G1xgreat*W*G1xless'...
      -G1yless*W*G1ygreat'+G1ygreat*W*G1yless'-G1zless*W*G1zgreat'+G1zgreat*W*G1zless');
    Ixx = -J^2/(2*pi^2).*(G1xless*W*G1xgreat'-G1xgreat*W*G1xless');
    Iyy = -J^2/(2*pi^2).*(G1yless*W*G1ygreat'-G1ygreat*W*G1yless');
    Izz = -J^2/(2*pi^2).*(G1zless*W*G1zgreat'-G1zgreat*W*G1zless');
    Ixz = -J^2/(4*pi^2).*(G1xless*W*G1zgreat'-G1xgreat*W*G1zless'+G1zless*W*G1xgreat'-G1zgreat*W*G1xless');
    Ixy = -J^2/(4*pi^2).*(G1xless*W*G1ygreat'-G1xgreat*W*G1yless'+G1yless*W*G1xgreat'-G1ygreat*W*G1xless');
    Iyz = -J^2/(4*pi^2).*(G1yless*W*G1zgreat'-G1ygreat*W*G1zless'+G1zless*W*G1ygreat'-G1zgreat*W*G1yless');
    Dx = -1i*J^2/(4*pi^2).*(G0less*W*G1xgreat'-G0great*W*G1xless'-G1xless*W*G0great'+G1xgreat*W*G0less');
    Dy = -1i*J^2/(4*pi^2).*(G0less*W*G1ygreat'-G0great*W*G1yless'-G1yless*W*G0great'+G1ygreat*W*G0less');
    Dz = -1i*J^2/(4*pi^2).*(G0less*W*G1zgreat'-G0great*W*G1zless'-G1zless*W*G0great'+G1zgreat*W*G0less');

    GJH = -J^2/(4*pi^2).*(G0less*WG*G0great'-G0great*WG*G0less'-G1xless*WG*G1xgreat'+G1xgreat*WG*G1xless'...
      -G1yless*WG*G1ygreat'+G1ygreat*WG*G1yless'-G1zless*WG*G1zgreat'+G1zgreat*WG*G1zless');
    GIxx = -J^2/(2*pi^2).*(G1xless*WG*G1xgreat'-G1xgreat*WG*G1xless');
    GIyy = -J^2/(2*pi^2).*(G1yless*WG*G1ygreat'-G1ygreat*WG*G1yless');
    GIzz = -J^2/(2*pi^2).*(G1zless*WG*G1zgreat'-G1zgreat*WG*G1zless');
    GIxz = -J^2/(4*pi^2).*(G1xless*WG*G1zgreat'-G1xgreat*WG*G1zless'+G1zless*WG*G1xgreat'-G1zgreat*WG*G1xless');
    GIxy = -J^2/(4*pi^2).*(G1xless*WG*G1ygreat'-G1xgreat*WG*G1yless'+G1yless*WG*G1xgreat'-G1ygreat*WG*G1xless');
    GIyz = -J^2/(4*pi^2).*(G1yless*WG*G1zgreat'-G1ygreat*WG*G1zless'+G1zless*WG*G1ygreat'-G1zgreat*WG*G1yless');
    GDx = -1i*J^2/(4*pi^2).*(G0less*WG*G1xgreat'-G0great*WG*G1xless'-G1xless*WG*G0great'+G1xgreat*WG*G0less');
    GDy = -1i*J^2/(4*pi^2).*(G0less*WG*G1ygreat'-G0great*WG*G1yless'-G1yless*WG*G0great'+G1ygreat*WG*G0less');
    GDz = -1i*J^2/(4*pi^2).*(G0less*WG*G1zgreat'-G0great*WG*G1zless'-G1zless*WG*G0great'+G1zgreat*WG*G0less');

    JH = real(JH)*dw^2;
    Ixx = real(Ixx)*dw^2;
    Iyy = real(Iyy)*dw^2;
    Izz = real(Izz)*dw^2;
    Ixz = real(Ixz)*dw^2;
    Ixy = real(Ixy)*dw^2;
    Iyz = real(Iyz)*dw^2;
    Dx = real(Dx)*dw^2;
    Dy = real(Dy)*dw^2;
    Dz = real(Dz)*dw^2;
    GJH = -imag(GJH)*dw^2;
    GIxx = -imag(GIxx)*dw^2;
    GIyy = -imag(GIyy)*dw^2;
    GIzz = -imag(GIzz)*dw^2;
    GIxz = -imag(GIxz)*dw^2;
    GIxy = -imag(GIxy)*dw^2;
    GIyz = -imag(GIyz)*dw^2;
    GDx = -imag(GDx)*dw^2;
    GDy = -imag(GDy)*dw^2;
    GDz = -imag(GDz)*dw^2;
end
