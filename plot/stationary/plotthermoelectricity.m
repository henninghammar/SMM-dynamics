xaxis = {epsvector};
xlabels = {'Gate voltage (1/Gamma)'};
labels = {'Conductance'; 'Spin conductance'; 'Charge Seebeck'; 'Spin Seebeck'; 'Heat conductance (electric)'; 'Figure of merit (charge)'; 'Figure of merit (spin)'; 'Peltier'};
plot_gate = [GL; GSL; ScL; SsL; kL; ZTL; ZTSL; T(1).*ScL];
plot_data = {plot_gate};
m = 0;
for parameter = 1:length(xlabels)
  for i = 1:length(labels)
    m = m + 1;
    h(m)=figure;
    plot_variables = plot_data{parameter};
    xaxis_variables = xaxis{parameter};
    plot(xaxis_variables, plot_variables(i,:))
    xlabel(xlabels(parameter))
    ylabel(labels(i))
  end
end
