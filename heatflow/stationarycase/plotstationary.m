
% xaxis = {eVvector; epsvector; Jvector; tempvector};
% xlabels = {'Bias voltage (1/Gamma)', 'Gate voltage (1/Gamma)', 'Exchange coupling (1/Gamma)', 'Temperature T_R (K) (T_L = 1 K)'};
% labels = {'Charge current', 'Spin current (z-direction)', 'Energy current', 'Heat current'};
% plot_bias = [Icbias; Iszbias; Iebias; Iqbias];
% plot_gate = [Icgate; Iszgate; Iegate; Iqgate];
% plot_J = [IcJ; IszJ; IeJ; IqJ];
% plot_temp = [Ictemp; Isztemp; Ietemp; Iqtemp];
% plot_data = {plot_bias; plot_gate; plot_J; plot_temp};

xaxis = {eVvector; epsvector};
xlabels = {'Bias voltage (1/Gamma)', 'Gate voltage (1/Gamma)'};
labels = {'Charge current', 'Spin current (z-direction)', 'Energy current', 'Heat current'};
plot_bias = [Icbias; Iszbias; Iebias; Iqbias];
plot_gate = [Icgate; Iszgate; Iegate; Iqgate];
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
% xaxis = {eVvector; epsvector; Jvector; tempvector};
% xlabels = {'Bias voltage (1/Gamma)', 'Gate voltage (1/Gamma)', 'Exchange coupling (1/Gamma)', 'Temperature T_R (K) (T_L = 1 K)'};
% Heisenberg = {JHbias; JHgate; JHJ; JHtemp};
% Ising = {Isingbias; Isinggate; IsingJ; Isingtemp};
% DM = {DMbias; DMgate; DMJ; DMtemp};
xaxis = {eVvector; epsvector};
xlabels = {'Bias voltage (1/Gamma)', 'Gate voltage (1/Gamma)'};
Heisenberg = {JHbias; JHgate};
Ising = {Isingbias; Isinggate};
DM = {DMbias; DMgate};
for parameter = 1:length(Heisenberg)
  i = i+1;
  h(i)=figure;
  hold on
  Heisenberg_values = Heisenberg{parameter};
  xaxis_variables = xaxis{parameter};
  plot(xaxis_variables, Heisenberg_values)
  xlabel(xlabels(parameter))
  ylabel('Heisenberg interactions')

  i = i+1;
  h(i)=figure;
  hold on
  Ising_values = Ising{parameter};
  for n = 1:6
    plot(xaxis_variables, Ising_values(n,:))
  end
  xlabel(xlabels(parameter))
  ylabel('Ising interactions')

  i = i+1;
  h(i)=figure;
  hold on
  DM_values = DM{parameter};
  for n = 1:3
    plot(xaxis_variables, DM_values(n,:))
  end
  xlabel(xlabels(parameter))
  ylabel('DM interactions')
end
