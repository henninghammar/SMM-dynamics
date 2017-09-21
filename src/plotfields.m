for loopen=1:length(values)
    pulse=values(loopen);
    inputFilename = sprintf(filename, outfolder, pulse);
    load(inputFilename)

    Bxmatrix(loopen,:)=real(SBx);
    mxmatrix(loopen,:)=real(Smx);
    ejxmatrix(loopen,:)=real(Sejx);
    jHxmatrix(loopen,:)=real(SjHx);
    jDMxmatrix(loopen,:)=real(SjDMx);
    jIxmatrix(loopen,:)=real(SjIx);
    Bymatrix(loopen,:)=real(SBy);
    mymatrix(loopen,:)=real(Smy);
    ejymatrix(loopen,:)=real(Sejy);
    jHymatrix(loopen,:)=real(SjHy);
    jDMymatrix(loopen,:)=real(SjDMy);
    jIymatrix(loopen,:)=real(SjIy);
    Bzmatrix(loopen,:)=real(SBz);
    mzmatrix(loopen,:)=real(Smz);
    ejzmatrix(loopen,:)=real(Sejz);
    jHzmatrix(loopen,:)=real(SjHz);
    jDMzmatrix(loopen,:)=real(SjDMz);
    jIzmatrix(loopen,:)=real(SjIz);
    Btotalx(loopen,:) = real(SjHx)+real(SjDMx)+real(SjIx)+real(SBx);
    Btotaly(loopen,:) = real(SjHy)+real(SjDMy)+real(SjIy)+real(SBy);
    Btotalz(loopen,:) = real(SjHz)+real(SjDMz)+real(SjIz)+real(SBz);

    jHmatrix(loopen,:)=real(jHt);

    Ixxmatrix(loopen,:)=real(jIxxt);
    Iyymatrix(loopen,:)=real(jIyyt);
    Izzmatrix(loopen,:)=real(jIzzt);

    DMxmatrix(loopen,:)=real(jDMxt);
    DMymatrix(loopen,:)=real(jDMyt);
    DMzmatrix(loopen,:)=real(jDMzt);
end

h30=figure(30);
contourf(t, values, Btotalx,100,'Linestyle','none')
title('Btotalx')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h30,strcat(outfolder,'beffx'),'fig');

h31=figure(31);
contourf(t, values, Btotaly,100,'Linestyle','none')
title('Btotaly')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h31,strcat(outfolder,'beffy'),'fig');

h32=figure(32);
contourf(t, values, Btotalz,100,'Linestyle','none')
title('Btotalz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h32,strcat(outfolder,'beffz'),'fig');

h1=figure(1);
contourf(t, values, Bxmatrix,100,'Linestyle','none')
title('Bx')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h1,strcat(outfolder,'Bxcontour'),'fig');

h33=figure(33);
contourf(t, values, ejxmatrix,100,'Linestyle','none')
title('ejx')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h33,strcat(outfolder,'ejxcontour'),'fig');

h34=figure(34);
contourf(t, values, mxmatrix,100,'Linestyle','none')
title('mx')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h34,strcat(outfolder,'mxcontour'),'fig');

h3=figure(3);
contourf(t, values, jHxmatrix,100,'Linestyle','none')
title('SjHx')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h3,strcat(outfolder,'SjHxcontour'),'fig');

h4=figure(4);
contourf(t, values, jIxmatrix,100,'Linestyle','none')
title('SIz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h4,strcat(outfolder,'SIxcontour'),'fig');

h5=figure(5);
contourf(t, values, jDMxmatrix,100,'Linestyle','none')
title('SDMx')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h5,strcat(outfolder,'SDMxcontour'),'fig');

h41=figure(41);
contourf(t, values, Bymatrix,100,'Linestyle','none')
title('By')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h41,strcat(outfolder,'Bycontour'),'fig');

h46=figure(46);
contourf(t, values, ejymatrix,100,'Linestyle','none')
title('ejy')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h46,strcat(outfolder,'ejycontour'),'fig');

h47=figure(47);
contourf(t, values, mymatrix,100,'Linestyle','none')
title('my')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h47,strcat(outfolder,'mycontour'),'fig');

h43=figure(43);
contourf(t, values, jHymatrix,100,'Linestyle','none')
title('SjHy')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h43,strcat(outfolder,'SjHycontour'),'fig');

h44=figure(44);
contourf(t, values, jIymatrix,100,'Linestyle','none')
title('SIy')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h44,strcat(outfolder,'SIycontour'),'fig');

h45=figure(45);
contourf(t, values, jDMymatrix,100,'Linestyle','none')
title('SDMy')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h45,strcat(outfolder,'SDMycontour'),'fig');

h51=figure(51);
contourf(t, values, Bzmatrix,100,'Linestyle','none')
title('Bz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h51,strcat(outfolder,'Bzcontour'),'fig');

h56=figure(56);
contourf(t, values, ejzmatrix,100,'Linestyle','none')
title('ejz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h56,strcat(outfolder,'ejzcontour'),'fig');

h57=figure(57);
contourf(t, values, mzmatrix,100,'Linestyle','none')
title('mz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h57,strcat(outfolder,'mzcontour'),'fig');

h53=figure(53);
contourf(t, values, jHzmatrix,100,'Linestyle','none')
title('SjHz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h53,strcat(outfolder,'SjHzcontour'),'fig');

h54=figure(54);
contourf(t, values, jIzmatrix,100,'Linestyle','none')
title('SIz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h54,strcat(outfolder,'SIzcontour'),'fig');

h55=figure(55);
contourf(t, values, jDMzmatrix,100,'Linestyle','none')
title('SDMz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h55,strcat(outfolder,'SDMzcontour'),'fig');
