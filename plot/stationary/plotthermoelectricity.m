
xaxis = {eVvector; epsvector};
xlabels = {'Bias voltage (1/Gamma)', 'Gate voltage (1/Gamma)'};
labels = {'Conductance'; 'Spin conductance'; 'Charge Seebeck'; 'Spin Seebeck'; 'Heat conductance (electric)'; 'Figure of merit (charge)'; 'Figure of merit (spin)'; 'Peltier'};
plot_gate = [GLgate; GSLgate; ScLgate; SsLgate; kLgate; ZTLgate; ZTSLgate; T(1).*ScLgate];
plot_bias = [GLbias; GSLbias; ScLbias; SsLbias; kLbias; ZTLbias; ZTSLbias; T(1).*ScLbias];
plot_data = {plot_bias; plot_gate};
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
%
% values = epsvector;
% labels = {'Transmission (up)', 'Transmission (down)'};
% m = 0;
% for i = 1:length(labels)
%   m = m + 1;
%   h(m)=figure;
%   hold on
%   for loop=1:length(values)
%     plot_data = {T_up(loop,:); T_down(loop,:)};
%     plot_variables = plot_data{i};
%     plot(plot_variables)
%     title(labels(i))
%     xlabel('Energy')
%     ylabel(labels(i))
%     legend()
%     %legend('eV = \Gamma', 'eV = 2\Gamma', 'eV = 3\Gamma', 'eV = 4\Gamma', 'eV = 5\Gamma');
%   end
% end
