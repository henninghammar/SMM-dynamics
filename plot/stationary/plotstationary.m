
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
labels = {'Charge current', 'Spin current (z-direction)', 'Energy current (L)', 'Heat current (L)', 'Particle current (L)', 'Energy current (R)', 'Heat current (R)', 'Particle current (R)', 'Net heat transfer', 'Net heat transfer (charge)', 'Net heat transfer (spin)'};
plot_gate = [Icgate; Iszgate; IeLgate; IqLgate; InLgate; IeRgate; IqRgate; InRgate; (IqLgate-IqRgate); (Iq0Lgate-Iq0Rgate); (Iq1Lgate-Iq1Rgate)];
plot_bias = [Icbias; Iszbias; IeLbias; IqLbias; InLbias; IeRbias; IqRbias; InRbias; (IqLbias-IqRbias); (Iq0Lbias-Iq0Rbias); (Iq1Lbias-Iq1Rbias)];
plot_data = {plot_bias; plot_gate};

% xaxis = {epsvector};
% xlabels = {'Gate voltage (1/Gamma)'};
% labels = {'Charge current', 'Spin current (z-direction)', 'Energy current (L)', 'Heat current (L)', 'Energy current (R)', 'Heat current (R)', 'Net heat transfer'};
% plot_gate = [Icgate; Iszgate; IeLgate; IqLgate; IeRgate; IqRgate; (IqLgate-IqRgate)];
% plot_data = {plot_gate};

m = 0;
for parameter = 1:length(xlabels)
  for i = 1:length(labels)
    m = m + 1;
    h(m)=figure(m);
    hold on
    plot_variables = plot_data{parameter};
    xaxis_variables = xaxis{parameter};
    plot(xaxis_variables, plot_variables(i,:))
    xlabel(xlabels(parameter))
    ylabel(labels(i))
  end
end
