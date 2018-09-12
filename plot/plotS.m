h1=figure(1);
hold on
plot(t, real(Sejx),'r')
plot(t, real(Smx))
plot(t, real(SjHx),'g')
plot(t, real(SjDMx),'black')
plot(t, real(SjIx))
title('Interactions x')
legend('ejx', 'mx', 'jHx', 'DMx', 'Ix')

h2=figure(2);
hold on
plot(t, real(Sejy),'r')
plot(t, real(Smy))
plot(t, real(SjHy),'g')
plot(t, real(SjDMy),'black')
plot(t, real(SjIy))
title('Interactions y')
legend('ejy', 'my', 'jHy', 'DMy', 'Iy')

h3=figure(3);
hold on
plot(t, real(Sejz),'r')
plot(t, real(Smz))
plot(t, real(SjHz),'g')
plot(t, real(SjDMz),'black')
plot(t, real(SjIz))
title('Interactions z')
legend('ejz', 'mz', 'jHz', 'DMz', 'Iz')

titles = {'ej', 'mv', 'Beff', 'jDM'};
plot_data = [ejvect; mvect; Beffvect; jDMvect];

titles = {'ej', 'Beff', 'jDM'};
plot_data = [ejvect; Beffvect; jDMvect];

for i = 1:3
h(i)=figure;
hold on
plot(real(plot_data(i+1,:)))
plot(real(plot_data(i+2,:)))
plot(real(plot_data(i+3,:)))
title(titles(i))
end

i = i+1;
h(i) = figure;
plot(real(jHt))
title('jH')

i = i+1;
h(i)=figure;
hold on
plot(real(jIvect(1,:)))
plot(real(jIvect(2,:)))
plot(real(jIvect(3,:)))
plot(real(jIvect(4,:)))
plot(real(jIvect(5,:)))
plot(real(jIvect(6,:)))
title('jI')
