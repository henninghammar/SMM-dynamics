h1=figure(4);
hold on
plot(real(SBx),'r')
plot(real(SjHx),'g')
plot(real(SjDMx),'black')
plot(real(SjIx))
title('Interactions x')
legend

h2=figure(5);
hold on
plot(real(SBy),'r')
plot(real(SjHy),'g')
plot(real(SjDMy),'black')
plot(real(SjIy))
title('Interactions y')
legend

h3=figure(6);
hold on
plot(real(SBz),'r')
plot(real(SjHz),'g')
plot(real(SjDMz),'black')
plot(real(SjIz))
title('Interactions z')
legend


% h2=figure(2);
% contourf(SjDMx)
% title('jDMx')
% colormap morgenstemning
% colorbar
%
% h3=figure(3);
% contourf(jDMytot)
% title('jDMy')
% colormap morgenstemning
% colorbar
%
% h4=figure(4);
% contourf(jDMztot)
% title('jDMz')
% colormap morgenstemning
% colorbar
%
% h5=figure(5);
% contourf(jIxxtot)
% title('jIxx')
% colormap morgenstemning
% colorbar
%
% h6=figure(6);
% contourf(jIxytot)
% title('jIxy')
% colormap morgenstemning
% colorbar
%
% h7=figure(7);
% contourf(jIxztot)
% title('jIxz')
% colormap morgenstemning
% colorbar
%
% h8=figure(8);
% contourf(jIyytot)
% title('jIyy')
% colormap morgenstemning
% colorbar
%
% h9=figure(9);
% contourf(jIyztot)
% title('jIyz')
% colormap morgenstemning
% colorbar
%
% h10=figure(10);
% contourf(jIzztot)
% title('jIzz')
% colormap morgenstemning
% colorbar
