h1=figure(1);
hold on
plot(real(SBx),'r')
plot(real(SjHx),'g')
plot(real(SjDMx),'black')
plot(real(SjIx))
title('Interactions x')
legend

h2=figure(2);
hold on
plot(real(SBy),'r')
plot(real(SjHy),'g')
plot(real(SjDMy),'black')
plot(real(SjIy))
title('Interactions y')
legend

h3=figure(3);
hold on
plot(real(SBz),'r')
plot(real(SjHz),'g')
plot(real(SjDMz),'black')
plot(real(SjIz))
title('Interactions z')
legend

%titles = {'ej', 'mv', 'Beff', 'jDM'};
%plot_data = [ejvect; mvect; Beffvect; jDMvect];
%
% titles = {'ej', 'Beff', 'jDM'};
% plot_data = [ejvect; Beffvect; jDMvect];
%
% for i = 1:3
% h(i)=figure;
% hold on
% plot(real(plot_data(i+1,:)))
% plot(real(plot_data(i+2,:)))
% plot(real(plot_data(i+3,:)))
% title(titles(i))
% end
%
% i = i+1;
% h(i) = figure;
% plot(real(jHt))
% title('jH')
%
% i = i+1;
% h(i)=figure;
% hold on
% plot(real(jIvect(1,:)))
% plot(real(jIvect(2,:)))
% plot(real(jIvect(3,:)))
% plot(real(jIvect(4,:)))
% plot(real(jIvect(5,:)))
% plot(real(jIvect(6,:)))
% title('jI')
