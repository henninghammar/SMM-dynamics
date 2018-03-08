for loop=1:length(values)
    gate1=values(loop);
    inputFilename = sprintf(filename, outfolder, gate1);
    load(inputFilename)

    Szm(loop,:)=Sz;
    ejzmatrix(loop,:)=real(Sejz);
end

% For gate
csteps = 100;
firststep = 0.5;
secondstep = 0.25;
thirdstep = 0.25;
%Light blue - green
mymap = [103*ones(firststep*csteps,1), 169*ones(firststep*csteps,1), 207*ones(firststep*csteps,1)
linspace(103,209,secondstep*csteps)', linspace(169,229,secondstep*csteps)', linspace(207,240,secondstep*csteps)'
linspace(209,255,thirdstep*csteps)', linspace(229,255,thirdstep*csteps)', linspace(240,255,thirdstep*csteps)'
linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;

% firststep = 0.5;
% secondstep = 0.5;
% thirdstep = 0.5;
% forthstep = 0.5;
% %Alt 2
% mymap = [255*ones(forthstep*csteps,1), 255*ones(forthstep*csteps,1), 255*ones(forthstep*csteps,1)
% linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
% linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
% zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;

h1=figure(1);
contourf(t, values, Szm,100,'Linestyle','none')
caxis([-1, 1])
title('Sz')
xlabel('time')
ylabel(axisvariable)
colormap(mymap)
colorbar
saveas(h1,strcat(outfolder,'Szcontour'),'fig');
%
%
% h2=figure(2);
% contourf(t, values, ejzmatrix,100,'Linestyle','none')
% caxis([-0.1, 0.1])
% title('ejz')
% xlabel('time')
% ylabel(axisvariable)
% colormap(mymap)
% colorbar
% saveas(h2,strcat(outfolder,'ejzcontour'),'fig');

% % For J
% csteps = 100;
% firststep = 0.5;
% secondstep = 0.3;
% thirdstep = 0.2;
%
% %Blue - green
% mymap1 = [33*ones(firststep*csteps,1), 102*ones(firststep*csteps,1), 171*ones(firststep*csteps,1)
% linspace(33,102,secondstep*csteps)', linspace(102,168,secondstep*csteps)', linspace(171,207,secondstep*csteps)'
% linspace(102,255,thirdstep*csteps)', linspace(168,255,thirdstep*csteps)', linspace(207,255,thirdstep*csteps)'
% linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
% linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
% zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;
%
% %Blue - red
% mymap2 = [33*ones(firststep*csteps,1), 102*ones(firststep*csteps,1), 171*ones(firststep*csteps,1)
% linspace(33,102,secondstep*csteps)', linspace(102,168,secondstep*csteps)', linspace(171,207,secondstep*csteps)'
% linspace(102,255,thirdstep*csteps)', linspace(168,255,thirdstep*csteps)', linspace(207,255,thirdstep*csteps)'
% linspace(255,244,thirdstep*csteps)', linspace(255,165,thirdstep*csteps)', linspace(255,130,thirdstep*csteps)'
% linspace(244,202,secondstep*csteps)', linspace(165,0,secondstep*csteps)', linspace(130,32,secondstep*csteps)'
% 202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;
%
% %Black - red
% mymap3 = [64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1)
% linspace(64,186,secondstep*csteps)', linspace(64,186,secondstep*csteps)', linspace(64,186,secondstep*csteps)'
% linspace(186,255,thirdstep*csteps)', linspace(186,255,thirdstep*csteps)', linspace(186,255,thirdstep*csteps)'
% linspace(255,244,thirdstep*csteps)', linspace(255,165,thirdstep*csteps)', linspace(255,130,thirdstep*csteps)'
% linspace(244,202,secondstep*csteps)', linspace(165,0,secondstep*csteps)', linspace(130,32,secondstep*csteps)'
% 202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;
%
% %Latte - blue
% mymap4 = [255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 204*ones(firststep*csteps,1)
% linspace(255,199,secondstep*csteps)', linspace(255,233,secondstep*csteps)', linspace(204,180,secondstep*csteps)'
% linspace(199,65,thirdstep*csteps)', linspace(233,182,thirdstep*csteps)', linspace(180,196,thirdstep*csteps)'
% linspace(65,34,thirdstep*csteps)', linspace(182,94,thirdstep*csteps)', linspace(196,168,thirdstep*csteps)'
% linspace(34,12,secondstep*csteps)', linspace(94,44,secondstep*csteps)', linspace(168,132,secondstep*csteps)'
% 12*ones(firststep*csteps,1), 44*ones(firststep*csteps,1), 132*ones(firststep*csteps,1)]./256;
%
% %Magenta - green
% mymap5 = [208*ones(firststep*csteps,1), 28*ones(firststep*csteps,1), 139*ones(firststep*csteps,1)
% linspace(208,241,secondstep*csteps)', linspace(28,182,secondstep*csteps)', linspace(139,218,secondstep*csteps)'
% linspace(241,247,thirdstep*csteps)', linspace(182,247,thirdstep*csteps)', linspace(218,247,thirdstep*csteps)'
% linspace(247,184,thirdstep*csteps)', linspace(247,225,thirdstep*csteps)', linspace(247,134,thirdstep*csteps)'
% linspace(184,77,secondstep*csteps)', linspace(225,172,secondstep*csteps)', linspace(134,38,secondstep*csteps)'
% 77*ones(firststep*csteps,1), 172*ones(firststep*csteps,1), 38*ones(firststep*csteps,1)]./256;
%
% %Blue - magenta
% mymap6 = [33*ones(firststep*csteps,1), 102*ones(firststep*csteps,1), 171*ones(firststep*csteps,1)
% linspace(33,102,secondstep*csteps)', linspace(102,168,secondstep*csteps)', linspace(171,207,secondstep*csteps)'
% linspace(102,255,thirdstep*csteps)', linspace(168,255,thirdstep*csteps)', linspace(207,255,thirdstep*csteps)'
% linspace(255,251,thirdstep*csteps)', linspace(255,180,thirdstep*csteps)', linspace(255,185,thirdstep*csteps)'
% linspace(251,247,secondstep*csteps)', linspace(180,104,secondstep*csteps)', linspace(185,161,secondstep*csteps)'
% 247*ones(firststep*csteps,1), 104*ones(firststep*csteps,1), 161*ones(firststep*csteps,1)]./256;
%
% %Pink - green
% mymap7 = [247*ones(firststep*csteps,1), 104*ones(firststep*csteps,1), 161*ones(firststep*csteps,1)
% linspace(247,251,secondstep*csteps)', linspace(104,180,secondstep*csteps)', linspace(161,185,secondstep*csteps)'
% linspace(251,255,thirdstep*csteps)', linspace(180,255,thirdstep*csteps)', linspace(185,255,thirdstep*csteps)'
% linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
% linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
% zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;
%
% %Light blue - green
% mymap8 = [103*ones(firststep*csteps,1), 169*ones(firststep*csteps,1), 207*ones(firststep*csteps,1)
% linspace(103,209,secondstep*csteps)', linspace(169,229,secondstep*csteps)', linspace(207,240,secondstep*csteps)'
% linspace(209,255,thirdstep*csteps)', linspace(229,255,thirdstep*csteps)', linspace(240,255,thirdstep*csteps)'
% linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
% linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
% zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;
%
% h1=figure(1);
% contourf(t, values, Szm,100,'Linestyle','none')
% hold on
% contour(t, values, Szm, [0.5 0.5] , 'linecolor', [0.2 0.2 0.2])
% contour(t, values, Szm, [-0.5 -0.5] , 'linecolor', [0.2 0.2 0.2]')
% contour(t, values, Szm, [0.25 0.25] , 'linecolor', [0.4 0.4 0.4]')
% contour(t, values, Szm, [-0.25 -0.25] , 'linecolor', [0.4 0.4 0.4]')
% contour(t, values, Szm, [0 0] , 'linecolor', [1 1 1])
% caxis([-1, 1])
% title('Sz')
% xlabel('time')
% ylabel(axisvariable)
% colormap(mymap8)
% colorbar
% saveas(h1,strcat(outfolder,'Szcontour'),'fig');
