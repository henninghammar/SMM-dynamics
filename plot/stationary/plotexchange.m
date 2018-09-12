xaxis = {eVvector; epsvector; Jvector; tempvector};
xlabels = {'Bias voltage (1/Gamma)', 'Gate voltage (1/Gamma)', 'Exchange coupling (1/Gamma)', 'Temperature T_R (K) (T_L = 1 K)'};
Heisenberg = {JHbias; JHgate; JHJ; JHtemp};
Ising = {Isingbias; Isinggate; IsingJ; Isingtemp};
DM = {DMbias; DMgate; DMJ; DMtemp};

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
