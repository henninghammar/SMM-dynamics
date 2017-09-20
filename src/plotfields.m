for i=1:length(values)
    pulse=values(i);
    inputFilename = sprintf(filename, outputFolder, pulse);
    load(inputFilename)

    ejxmatrix(i,:)=real(Sejx);
    jHxmatrix(i,:)=real(SjHx);
    jDMxmatrix(i,:)=real(SjDMx);
    jIxmatrix(i,:)=real(SjIx);
    ejymatrix(i,:)=real(Sejy);
    jHymatrix(i,:)=real(SjHy);
    jDMymatrix(i,:)=real(SjDMy);
    jIymatrix(i,:)=real(SjIy);
    ejzmatrix(i,:)=real(Sejz);
    jHzmatrix(i,:)=real(SjHz);
    jDMzmatrix(i,:)=real(SjDMz);
    jIzmatrix(i,:)=real(SjIz);
    Beffectivex(i,:) = real(SjHx)+real(SjDMx)+real(SjIx)+real(SBx);
    Beffectivey(i,:) = real(SjHy)+real(SjDMy)+real(SjIy)+real(SBy);
    Beffectivez(i,:) = real(SjHz)+real(SjDMz)+real(SjIz)+real(SBz);

    jHmatrix(i,:)=real(jHt);

    Ixxmatrix(i,:)=real(jIxxt);
    Iyymatrix(i,:)=real(jIyyt);
    Izzmatrix(i,:)=real(jIzzt);

    DMxmatrix(i,:)=real(jDMxt);
    DMymatrix(i,:)=real(jDMyt);
    DMzmatrix(i,:)=real(jDMzt);
end
%

h30=figure(30);
contourf(t, values, Beffectivex,100,'Linestyle','none')
title('Beffx')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h30,'beffx','fig');

h31=figure(31);
contourf(t, values, Beffectivey,100,'Linestyle','none')
title('Beffy')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h31,'beffy','fig');

h32=figure(32);
contourf(t, values, Beffectivez,100,'Linestyle','none')
title('Beffz')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h32,'beffz','fig');

h1=figure(1);
contourf(t, values, ejxmatrix,100,'Linestyle','none')
title('ejx')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h1,'ejxcontour','fig');

h3=figure(3);
contourf(t, values, jHxmatrix,100,'Linestyle','none')
title('SjHx')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h3,'SjHxcontour','fig');

h4=figure(4);
contourf(t, values, jIxmatrix,100,'Linestyle','none')
title('SIz')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h4,'SIxcontour','fig');

h5=figure(5);
contourf(t, values, jDMxmatrix,100,'Linestyle','none')
title('SDMx')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h5,'SDMxcontour','fig');

h41=figure(41);
contourf(t, values, ejymatrix,100,'Linestyle','none')
title('ejy')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h41,'ejycontour','fig');

h43=figure(43);
contourf(t, values, jHymatrix,100,'Linestyle','none')
title('SjHy')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h43,'SjHycontour','fig');

h44=figure(44);
contourf(t, values, jIymatrix,100,'Linestyle','none')
title('SIy')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h44,'SIycontour','fig');

h45=figure(45);
contourf(t, values, jDMymatrix,100,'Linestyle','none')
title('SDMy')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h45,'SDMycontour','fig');

h51=figure(51);
contourf(t, values, ejzmatrix,100,'Linestyle','none')
title('ejz')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h51,'ejzcontour','fig');

h53=figure(53);
contourf(t, values, jHzmatrix,100,'Linestyle','none')
title('SjHz')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h53,'SjHzcontour','fig');

h54=figure(54);
contourf(t, values, jIzmatrix,100,'Linestyle','none')
title('SIz')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h54,'SIzcontour','fig');

h55=figure(55);
contourf(t, values, jDMzmatrix,100,'Linestyle','none')
title('SDMz')
xlabel('T(s)')
ylabel('Pulse length')
colormap morgenstemning
colorbar
saveas(h55,'SDMzcontour','fig');
