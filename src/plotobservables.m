for loop=1:length(values)
    gate1=values(loop);
    inputFilename = sprintf(filename, outfolder, gate1);
    load(inputFilename)

    Szn(loop)=Sz(end);
    Icn(loop)=sum(Ic)/length(Ic1);
    ISzn(loop)=sum(Isz)/length(Isz0);

    Sxm(loop,:)=Sx;
    Sym(loop,:)=Sy;
    Szm(loop,:)=Sz;
    Icm(loop,:)=Ic;
    ISxm(loop,:)=Isx;
    ISym(loop,:)=Isy;
    ISzm(loop,:)=Isz;
end

Sxm=round(Sxm.*10^3)./10^3;
Sym=round(Sym.*10^3)./10^3;
Szm=round(Szm.*10^3)./10^3;

h10=figure(10);
plot(values,Icn)
title('Charge current')
xlabel(axisvariable)
ylabel('Ic')
saveas(h10,strcat(outfolder,'Ic'),'fig');

h11=figure(11);
plot(values,ISzn)
title('Spin current')
xlabel(axisvariable)
ylabel('Isz')
saveas(h11,strcat(outfolder,'Isz'),'fig');

h12=figure(12);
plot(values,Szn)
title('Spin Sz-direction')
xlabel(axisvariable)
ylabel('Sz')
saveas(h12,strcat(outfolder,'Sz'),'fig');


h13=figure(13);
contourf(t, values, Sxm,100,'Linestyle','none')
title('Sx')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h13,strcat(outfolder,'Sxcontour'),'fig');

h14=figure(14);
contourf(t, values, Sym,100,'Linestyle','none')
title('Sy')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h14,strcat(outfolder,'Sycontour'),'fig');

h15=figure(15);
contourf(t, values, Szm,100,'Linestyle','none')
title('Sz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h15,strcat(outfolder,'Szcontour'),'fig');

h16=figure(16);
contourf(t, values, Icm,100,'Linestyle','none')
title('Ic')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h16,strcat(outfolder,'Iccontour'),'fig');

h17=figure(17);
contourf(t, values, ISxm,100,'Linestyle','none')
title('Isx')
xlabel('time')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h17,strcat(outfolder,'Isxcontour'),'fig');


h18=figure(18);
contourf(t, values, ISym,100,'Linestyle','none')
title('Isy')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h18,strcat(outfolder,'Isycontour'),'fig');


h19=figure(19);
contourf(t, values, ISzm,100,'Linestyle','none')
title('Isz')
xlabel('time')
ylabel(axisvariable)
colormap morgenstemning
colorbar
saveas(h19,strcat(outfolder,'Iszcontour'),'fig');
