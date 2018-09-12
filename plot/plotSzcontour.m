for loop=1:length(values)
    gate1=values(loop);
    inputFilename = sprintf(filename, outfolder, gate1);
    load(inputFilename)

    Szm(loop,:)=Sz;
    ejzmatrix(loop,:)=real(Sejz);
    Iszm(loop,:)=real(Isz);
end

mu=eV(1)-eV(2);
%values = mu*values/(2*pi);

csteps = 100;

%Different levels
firststep = 0.3;
secondstep = 0.7;
thirdstep = 0.5;

% Different color schemes
%Blue - green
mymap1 = [33*ones(firststep*csteps,1), 102*ones(firststep*csteps,1), 171*ones(firststep*csteps,1)
linspace(33,102,secondstep*csteps)', linspace(102,168,secondstep*csteps)', linspace(171,207,secondstep*csteps)'
linspace(102,255,thirdstep*csteps)', linspace(168,255,thirdstep*csteps)', linspace(207,255,thirdstep*csteps)'
linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;
%
%Blue - red
mymap2 = [33*ones(firststep*csteps,1), 102*ones(firststep*csteps,1), 171*ones(firststep*csteps,1)
linspace(33,102,secondstep*csteps)', linspace(102,168,secondstep*csteps)', linspace(171,207,secondstep*csteps)'
linspace(102,255,thirdstep*csteps)', linspace(168,255,thirdstep*csteps)', linspace(207,255,thirdstep*csteps)'
linspace(255,244,thirdstep*csteps)', linspace(255,165,thirdstep*csteps)', linspace(255,130,thirdstep*csteps)'
linspace(244,202,secondstep*csteps)', linspace(165,0,secondstep*csteps)', linspace(130,32,secondstep*csteps)'
202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;
%
% %Black - red
mymap3 = [64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1)
linspace(64,186,secondstep*csteps)', linspace(64,186,secondstep*csteps)', linspace(64,186,secondstep*csteps)'
linspace(186,255,thirdstep*csteps)', linspace(186,255,thirdstep*csteps)', linspace(186,255,thirdstep*csteps)'
linspace(255,244,thirdstep*csteps)', linspace(255,165,thirdstep*csteps)', linspace(255,130,thirdstep*csteps)'
linspace(244,202,secondstep*csteps)', linspace(165,0,secondstep*csteps)', linspace(130,32,secondstep*csteps)'
202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;

%Latte - blue
mymap4 = [255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 204*ones(firststep*csteps,1)
linspace(255,199,secondstep*csteps)', linspace(255,233,secondstep*csteps)', linspace(204,180,secondstep*csteps)'
linspace(199,65,thirdstep*csteps)', linspace(233,182,thirdstep*csteps)', linspace(180,196,thirdstep*csteps)'
linspace(65,34,thirdstep*csteps)', linspace(182,94,thirdstep*csteps)', linspace(196,168,thirdstep*csteps)'
linspace(34,12,secondstep*csteps)', linspace(94,44,secondstep*csteps)', linspace(168,132,secondstep*csteps)'
12*ones(firststep*csteps,1), 44*ones(firststep*csteps,1), 132*ones(firststep*csteps,1)]./256;

%Magenta - green
mymap5 = [208*ones(firststep*csteps,1), 28*ones(firststep*csteps,1), 139*ones(firststep*csteps,1)
linspace(208,241,secondstep*csteps)', linspace(28,182,secondstep*csteps)', linspace(139,218,secondstep*csteps)'
linspace(241,247,thirdstep*csteps)', linspace(182,247,thirdstep*csteps)', linspace(218,247,thirdstep*csteps)'
linspace(247,184,thirdstep*csteps)', linspace(247,225,thirdstep*csteps)', linspace(247,134,thirdstep*csteps)'
linspace(184,77,secondstep*csteps)', linspace(225,172,secondstep*csteps)', linspace(134,38,secondstep*csteps)'
77*ones(firststep*csteps,1), 172*ones(firststep*csteps,1), 38*ones(firststep*csteps,1)]./256;

%Blue - magenta
mymap6 = [33*ones(firststep*csteps,1), 102*ones(firststep*csteps,1), 171*ones(firststep*csteps,1)
linspace(33,102,secondstep*csteps)', linspace(102,168,secondstep*csteps)', linspace(171,207,secondstep*csteps)'
linspace(102,255,thirdstep*csteps)', linspace(168,255,thirdstep*csteps)', linspace(207,255,thirdstep*csteps)'
linspace(255,251,thirdstep*csteps)', linspace(255,180,thirdstep*csteps)', linspace(255,185,thirdstep*csteps)'
linspace(251,247,secondstep*csteps)', linspace(180,104,secondstep*csteps)', linspace(185,161,secondstep*csteps)'
247*ones(firststep*csteps,1), 104*ones(firststep*csteps,1), 161*ones(firststep*csteps,1)]./256;

%Pink - green
mymap7 = [247*ones(firststep*csteps,1), 104*ones(firststep*csteps,1), 161*ones(firststep*csteps,1)
linspace(247,251,secondstep*csteps)', linspace(104,180,secondstep*csteps)', linspace(161,185,secondstep*csteps)'
linspace(251,255,thirdstep*csteps)', linspace(180,255,thirdstep*csteps)', linspace(185,255,thirdstep*csteps)'
linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;

%Light blue - green
mymap8 = [103*ones(firststep*csteps,1), 169*ones(firststep*csteps,1), 207*ones(firststep*csteps,1)
linspace(103,209,secondstep*csteps)', linspace(169,229,secondstep*csteps)', linspace(207,240,secondstep*csteps)'
linspace(209,255,thirdstep*csteps)', linspace(229,255,thirdstep*csteps)', linspace(240,255,thirdstep*csteps)'
linspace(255,89,thirdstep*csteps)', linspace(255,181,thirdstep*csteps)', linspace(255,171,thirdstep*csteps)'
linspace(89,0,secondstep*csteps)', linspace(181,102,secondstep*csteps)', linspace(171,94,secondstep*csteps)'
zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;

%Light blue to white to green
mymap9 = [103*ones(firststep*csteps,1), 169*ones(firststep*csteps,1), 207*ones(firststep*csteps,1)
linspace(103,255,secondstep*csteps)', linspace(169,255,secondstep*csteps)', linspace(207,255,secondstep*csteps)'
255*ones(thirdstep*csteps,1), 255*ones(thirdstep*csteps,1), 255*ones(thirdstep*csteps,1)
linspace(255,0,secondstep*csteps)', linspace(255,102,secondstep*csteps)', linspace(255,94,secondstep*csteps)'
zeros(firststep*csteps,1), 102*ones(firststep*csteps,1), 94*ones(firststep*csteps,1)]./256;

% %Grey - red
mymap10 = [186*ones(firststep*csteps,1), 186*ones(firststep*csteps,1), 186*ones(firststep*csteps,1)
linspace(186,255,(secondstep+thirdstep)*csteps)', linspace(186,255,(secondstep+thirdstep)*csteps)', linspace(186,255,(secondstep+thirdstep)*csteps)'
linspace(255,244,thirdstep*csteps)', linspace(255,165,thirdstep*csteps)', linspace(255,130,thirdstep*csteps)'
linspace(244,202,secondstep*csteps)', linspace(165,0,secondstep*csteps)', linspace(130,32,secondstep*csteps)'
202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;

% %White - red
mymap11 = [255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1)
linspace(255,202,secondstep*csteps)', linspace(255,0,secondstep*csteps)', linspace(255,32,secondstep*csteps)'
202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;

%Blue red whithout white in between
mymap12 = [33*ones(firststep*csteps,1), 102*ones(firststep*csteps,1), 171*ones(firststep*csteps,1)
linspace(33,102,secondstep*csteps)', linspace(102,168,secondstep*csteps)', linspace(171,207,secondstep*csteps)'
linspace(102,244,thirdstep*csteps)', linspace(168,165,thirdstep*csteps)', linspace(207,130,thirdstep*csteps)'
linspace(244,202,secondstep*csteps)', linspace(165,0,secondstep*csteps)', linspace(130,32,secondstep*csteps)'
202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;

%Light blue to white to red
mymap13 = [103*ones(firststep*csteps,1), 169*ones(firststep*csteps,1), 207*ones(firststep*csteps,1)
linspace(103,255,secondstep*csteps)', linspace(169,255,secondstep*csteps)', linspace(207,255,secondstep*csteps)'
255*ones(thirdstep*csteps,1), 255*ones(thirdstep*csteps,1), 255*ones(thirdstep*csteps,1)
linspace(255,202,secondstep*csteps)', linspace(255,0,secondstep*csteps)', linspace(255,32,secondstep*csteps)'
202*ones(firststep*csteps,1), zeros(firststep*csteps,1), 32*ones(firststep*csteps,1)]./256;

% %Black to red to white
mymap14 = [64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1)
linspace(64,202,secondstep*csteps)', linspace(64,0,secondstep*csteps)', linspace(64,32,secondstep*csteps)'
linspace(202,255,secondstep*csteps)', linspace(0,255,secondstep*csteps)', linspace(32,255,secondstep*csteps)'
255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1)]./256;

% %Black to latte to white
mymap15 = [64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1), 64*ones(firststep*csteps,1)
linspace(64,244,secondstep*csteps)', linspace(64,165,secondstep*csteps)', linspace(64,130,secondstep*csteps)'
linspace(244,255,secondstep*csteps)', linspace(165,255,secondstep*csteps)', linspace(130,255,secondstep*csteps)'
255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1)]./256;

% %Black to pink to white
mymap16 = [37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1)
linspace(37,223,secondstep*csteps)', linspace(37,101,secondstep*csteps)', linspace(37,176,secondstep*csteps)'
linspace(223,255,secondstep*csteps)', linspace(101,255,secondstep*csteps)', linspace(176,255,secondstep*csteps)'
255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1)]./256;

% %Black to bright pink to white
mymap17 = [37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1)
linspace(37,221,secondstep*csteps)', linspace(37,28,secondstep*csteps)', linspace(37,119,secondstep*csteps)'
linspace(221,255,secondstep*csteps)', linspace(28,255,secondstep*csteps)', linspace(119,255,secondstep*csteps)'
255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1)]./256;

% %Black to yellow to white
mymap18 = [37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1)
linspace(37,254,secondstep*csteps)', linspace(37,204,secondstep*csteps)', linspace(37,41,secondstep*csteps)'
linspace(254,255,secondstep*csteps)', linspace(204,255,secondstep*csteps)', linspace(41,255,secondstep*csteps)'
255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1)]./256;

% %Black to orange to white
mymap19 = [37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1), 37*ones(firststep*csteps,1)
linspace(37,253,secondstep*csteps)', linspace(37,141,secondstep*csteps)', linspace(37,60,secondstep*csteps)'
linspace(253,255,secondstep*csteps)', linspace(141,255,secondstep*csteps)', linspace(60,255,secondstep*csteps)'
255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1), 255*ones(firststep*csteps,1)]./256;

h1=figure(1);
contourf(t, values, Szm,100,'Linestyle','none')
caxis([-1, 1])
title('Sz')
xlabel('time')
ylabel(axisvariable)
colormap(mymap18)
colorbar
#saveas(h1,strcat(outfolder,'Szcontour'),'fig');

h2=figure(2);
contourf(t, values, ejzmatrix,100,'Linestyle','none')
caxis([-0.1, 0.1])
title('ejz')
xlabel('time')
ylabel(axisvariable)
colormap(mymap)
colorbar
% saveas(h2,strcat(outfolder,'ejzcontour'),'fig');

h3=figure(3);
contourf(t, values, Iszm,100,'Linestyle','none')
%caxis([-1, 1])
title('Spin current (z-direction)')
xlabel('time')
ylabel(axisvariable)
colormap(mymap18)
colorbar
