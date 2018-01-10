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
    Btotalx(loopen,:) = real(SjHx)+real(SjDMx)+real(SjIx)+real(SBx)+real(Smx);
    Btotaly(loopen,:) = real(SjHy)+real(SjDMy)+real(SjIy)+real(SBy)+real(Smy);
    Btotalz(loopen,:) = real(SjHz)+real(SjDMz)+real(SjIz)+real(SBz)+real(Smz);

    jHmatrix(loopen,:)=real(jHt);

    Ixxmatrix(loopen,:)=real(jIvect(1,:));
    Iyymatrix(loopen,:)=real(jIvect(2,:));
    Izzmatrix(loopen,:)=real(jIvect(3,:));

    DMxmatrix(loopen,:)=real(jDMvect(1,:));
    DMymatrix(loopen,:)=real(jDMvect(2,:));
    DMzmatrix(loopen,:)=real(jDMvect(3,:));
end
% %
% % h1=figure(1);
% % contourf(t, values, Btotalx,100,'Linestyle','none')
% % title('Btotalx')
% % xlabel('time')
% % ylabel(axisvariable)
% % colormap morgenstemning
% % colorbar
% % saveas(h1,strcat(outfolder,'beffx'),'fig');
% %
% % h2=figure(2);
% % contourf(t, values, Btotaly,100,'Linestyle','none')
% % title('Btotaly')
% % xlabel('time')
% % ylabel(axisvariable)
% % colormap morgenstemning
% % colorbar
% % saveas(h2,strcat(outfolder,'beffy'),'fig');
% %
% % h3=figure(3);
% % contourf(t, values, Btotalz,100,'Linestyle','none')
% % title('Btotalz')
% % xlabel('time')
% % ylabel(axisvariable)
% % colormap morgenstemning
% % colorbar
% % saveas(h3,strcat(outfolder,'beffz'),'fig');
%
% h11=figure(11);
% contourf(t, values, Bxmatrix,100,'Linestyle','none')
% title('Bx')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h11,strcat(outfolder,'Bxcontour'),'fig');
%
% h12=figure(12);
% contourf(t, values, ejxmatrix,100,'Linestyle','none')
% title('ejx')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h12,strcat(outfolder,'ejxcontour'),'fig');
%
% h13=figure(13);
% contourf(t, values, mxmatrix,100,'Linestyle','none')
% title('mx')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h13,strcat(outfolder,'mxcontour'),'fig');
%
% h14=figure(14);
% contourf(t, values, jHxmatrix,100,'Linestyle','none')
% title('SjHx')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h14,strcat(outfolder,'SjHxcontour'),'fig');
%
% h15=figure(15);
% contourf(t, values, jIxmatrix,100,'Linestyle','none')
% title('SIx')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h15,strcat(outfolder,'SIxcontour'),'fig');
%
% h16=figure(16);
% contourf(t, values, jDMxmatrix,100,'Linestyle','none')
% title('SDMx')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h16,strcat(outfolder,'SDMxcontour'),'fig');
% %
% % h17=figure(17);
% % contourf(t, values, DMxmatrix,100,'Linestyle','none')
% % title('DMx')
% % xlabel('time')
% % ylabel(axisvariable)
% % colormap morgenstemning
% % colorbar
% % saveas(h17,strcat(outfolder,'DMxcontour'),'fig');
% %
% % h18=figure(18);
% % contourf(t, values, Ixxmatrix,100,'Linestyle','none')
% % title('Ixx')
% % xlabel('time')
% % ylabel(axisvariable)
% % colormap morgenstemning
% % colorbar
% % saveas(h18,strcat(outfolder,'Ixxcontour'),'fig');
%
% h21=figure(21);
% contourf(t, values, Bymatrix,100,'Linestyle','none')
% title('By')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h21,strcat(outfolder,'Bycontour'),'fig');
%
% h22=figure(22);
% contourf(t, values, ejymatrix,100,'Linestyle','none')
% title('ejy')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h22,strcat(outfolder,'ejycontour'),'fig');
%
% h23=figure(23);
% contourf(t, values, mymatrix,100,'Linestyle','none')
% title('my')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h23,strcat(outfolder,'mycontour'),'fig');
%
% h24=figure(24);
% contourf(t, values, jHymatrix,100,'Linestyle','none')
% title('SjHy')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h24,strcat(outfolder,'SjHycontour'),'fig');
%
% h25=figure(25);
% contourf(t, values, jIymatrix,100,'Linestyle','none')
% title('SIy')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h25,strcat(outfolder,'SIycontour'),'fig');
%
% h26=figure(26);
% contourf(t, values, jDMymatrix,100,'Linestyle','none')
% title('SDMy')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h26,strcat(outfolder,'SDMycontour'),'fig');
% %
% % h27=figure(27);
% % contourf(t, values, DMymatrix,100,'Linestyle','none')
% % title('DMy')
% % xlabel('time')
% % ylabel(axisvariable)
% % colormap morgenstemning
% % colorbar
% % saveas(h27,strcat(outfolder,'DMycontour'),'fig');
% %
% % h28=figure(28);
% % contourf(t, values, Iyymatrix,100,'Linestyle','none')
% % title('Iyy')
% % xlabel('time')
% % ylabel(axisvariable)
% % colormap morgenstemning
% % colorbar
% % saveas(h28,strcat(outfolder,'Iyycontour'),'fig');
%
% h31=figure(31);
% contourf(t, values, Bzmatrix,100,'Linestyle','none')
% title('Bz')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h31,strcat(outfolder,'Bzcontour'),'fig');

h32=figure(32);
contourf(t, values, ejzmatrix,100,'Linestyle','none')
title('ejz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h32,strcat(outfolder,'ejzcontour'),'fig');

h33=figure(33);
contourf(t, values, mzmatrix,100,'Linestyle','none')
title('mz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h33,strcat(outfolder,'mzcontour'),'fig');

h34=figure(34);
contourf(t, values, jHzmatrix,100,'Linestyle','none')
title('SjHz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h34,strcat(outfolder,'SjHzcontour'),'fig');

h35=figure(35);
contourf(t, values, jIzmatrix,100,'Linestyle','none')
title('SIz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h35,strcat(outfolder,'SIzcontour'),'fig');

h36=figure(36);
contourf(t, values, jDMzmatrix,100,'Linestyle','none')
title('SDMz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h36,strcat(outfolder,'SDMzcontour'),'fig');
%
% h37=figure(37);
% contourf(t, values, DMzmatrix,100,'Linestyle','none')
% title('DMz')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h37,strcat(outfolder,'DMzcontour'),'fig');
%
% h38=figure(38);
% contourf(t, values, Izzmatrix,100,'Linestyle','none')
% title('Izz')
% xlabel('time')
% ylabel(axisvariable)
% colormap morgenstemning
% colorbar
% saveas(h38,strcat(outfolder,'Izzcontour'),'fig');
