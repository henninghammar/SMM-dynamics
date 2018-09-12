filename = 'test';

labels = {'Spin (z-direction)', 'Heat current (spin) (L)', 'Heat current (spin) (R)'};
plot_data = {Sz; IqL1; IqR1};
m = 1;
h(m) = figure;
for i = 1:length(labels)
  hold on
  plot_variables = plot_data{i};
  plot(t, plot_variables)
  xlabel('Time')
end
legend(labels)
saveas(h(m),strcat(filename,'heatcurrent'),'fig');
m = m +1;

labels = {'Spin (z-direction)', 'Heat current (L)', 'Heat current (R)'};
plot_data = {Sz; IqL; IqR};
h(m) = figure;
for i = 1:length(labels)
  hold on
  plot_variables = plot_data{i};
  plot(t, plot_variables)
  xlabel('Time')
end
legend(labels)
saveas(h(m),strcat(filename,'heatcurrent'),'fig');
m = m +1;

labels = {'Net heat flow (charge)', 'Net heat flow (spin)', 'Net heat flow (total)'};
plot_data = {(IqL0-IqR0); (IqL1-IqR1); (IqL-IqR)};
h(m) = figure;
for i = 1:length(labels)
  hold on
  plot_variables = plot_data{i};
  plot(t, plot_variables)
  xlabel('Time')
end
legend(labels)
saveas(h(m),strcat(filename,'heatflow'),'fig');

labels = {'Net particle flow (charge)', 'Net particle flow (spin)', 'Net particle flow (total)'};
plot_data = {(InL0+InR0); (InL1+InR1); (InL+InR)};
h(m) = figure;
for i = 1:length(labels)
  hold on
  plot_variables = plot_data{i};
  plot(t, plot_variables)
  xlabel('Time')
end
legend(labels)
saveas(h(m),strcat(filename,'netparticleflow'),'fig');

nup = (n/2 + mvect(3,:))/2;
ndown = (n/2 - mvect(3,:))/2;

mentropyTot = -nup.*log(nup)-ndown.*log(ndown);
mentropyDiff = -nup.*log(nup)+ndown.*log(ndown);

labels = {'Spin up-occupation', 'Spin down-occupation', 'Magnetic entropy'};
plot_data = {nup; ndown; mentropyTot};
h(m) = figure;
for i = 1:length(labels)
  hold on
  plot_variables = plot_data{i};
  plot(t, plot_variables)
  xlabel('Time')
end
legend(labels)
saveas(h(m),strcat(filename,'entropy'),'fig');
