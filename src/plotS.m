% h1=figure(1);
% hold on
% plot(real(SBx),'r')
% plot(real(SjHx),'g')
% plot(real(SjDMx),'black')
% plot(real(SjIx))
% title('Interactions x')
% legend
%
% h2=figure(2);
% hold on
% plot(real(SBy),'r')
% plot(real(SjHy),'g')
% plot(real(SjDMy),'black')
% plot(real(SjIy))
% title('Interactions y')
% legend
%
% h3=figure(3);
% hold on
% plot(real(SBz),'r')
% plot(real(SjHz),'g')
% plot(real(SjDMz),'black')
% plot(real(SjIz))
% title('Interactions z')
% legend

%jHt(j)=trapz(tau,jH);
%jIvect(:,j)=[trapz(tau,jIxx), trapz(tau,jIyy), trapz(tau,jIzz), trapz(tau,jIxy), trapz(tau,jIxz), trapz(tau,jIyz)]

titles = ['ej', 'mv', 'Beff', 'jDM'];
plot_data = [ejvect; mvect; Beffvect; jDMvect];

for i = 1:3
h(i)=figure(i);
hold on
plot(real(plot_data(i+1,:)))
plot(real(plot_data(i+2,:)))
plot(real(plot_data(i+3,:)))
title(titles(i))
legend
end
