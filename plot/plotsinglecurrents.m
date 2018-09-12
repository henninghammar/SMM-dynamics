labels = {'Charge current', 'Spin current (z-direction)', 'Energy current', 'Heat current'};
plot_data = {Ic; Isz; IeL; IqL};
m = 0;
for i = 1:length(labels)
  m = m + 1;
  h(m)=figure;
  plot_variables = plot_data{i};
  plot(t, plot_variables)
  xlabel('Time')
  ylabel(labels(i))
end
