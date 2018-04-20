for loop=1:length(values)
    gate1=values(loop);
    inputFilename = sprintf(filename, outfolder, gate1);
    load(inputFilename)

    Sxm(loop,:)=Sx;
    Sym(loop,:)=Sy;
    Szm(loop,:)=Sz;
    Icm(loop,:)=Ic;
    ISxm(loop,:)=Isx;
    ISym(loop,:)=Isy;
    ISzm(loop,:)=Isz;
    IeLm(loop,:)=IeL1;
    IqLm(loop,:)=IqL1;
end

labels = {'Spin (z-direction)', 'Charge current', 'Spin current (z-direction)', 'Energy current', 'Heat current'};
names = {'Szcontour', 'Iccontour', 'Iszcontour', 'Iecontour', 'Iqcontour'};
plot_data = {Szm; Icm; ISzm; IeLm; IqLm};
m = 0;
for i = 1:length(labels)
  m = m + 1;
  h(m)=figure;
  hold on
  plot_variables = plot_data{i};
  contourf(t, values, plot_variables, 100,'Linestyle','none')
  title(labels(i))
  xlabel('Time')
  ylabel(axisvariable)
  colormap morgenstemning
  colorbar
  saveas(h(m),strcat(outfolder,names{i}),'fig');
end

labels = {'Spin (z-direction)', 'Charge current', 'Spin current (z-direction)', 'Energy current', 'Heat current'};
names = {'Szend', 'Icend', 'Iszend', 'Ieend', 'Iqend'};
plot_data = {Szm(:,end); Icm(:,end); ISzm(:,end); IeLm(:,end); IqLm(:,end)};
m = 0;
for i = 1:length(labels)
  m = m + 1;
  h(m)=figure;
  hold on
  plot_variables = plot_data{i};
  plot(values, plot_variables)
  title(labels(i))
  xlabel(axisvariable)
  ylabel(labels(i))
  saveas(h(m),strcat(outfolder,names{i}),'fig');
end
