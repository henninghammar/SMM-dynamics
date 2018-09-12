#Add in your simulation folder
clear

addpath('../../../plot')

values=[1:5];
filename = '%s/...%02d.mat';
outfolder = 'output/.../';
axisvariable = 'Bias voltage';

plotcurrents

plotobservables

plotSzcontour

plotfields
