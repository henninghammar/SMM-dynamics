h1=figure(1);
hold on
plot(t,real(SBx),'r')
plot(t,real(SjHx),'g')
plot(t,real(SjDMx),'black')
plot(t,real(SjIx))
plot(t,real(SBx)+real(SjHx)+real(SjDMx)+real(SjIx),'black','linewidth', 2)
title('Interactions x')
legend

h2=figure(2);
hold on
plot(t,real(SBy),'r')
plot(t,real(SjHy),'g')
plot(t,real(SjDMy),'black')
plot(t,real(SjIy))
plot(t,real(SBy)+real(SjHy)+real(SjDMy)+real(SjIy),'black','linewidth', 2)
title('Interactions y')
legend

h3=figure(3);
hold on
ylim([-0.4,0.4])
plot(t,real(SBz),'r')
plot(t,real(SjHz),'g')
plot(t,real(SjDMz),'black')
plot(t,real(SjIz))
plot(t,real(SBz)+real(SjHz)+real(SjDMz)+real(SjIz),'black','linewidth', 2)
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
