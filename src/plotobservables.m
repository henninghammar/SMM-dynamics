for j=1:length(value)
    gate1=value(j);
    inputFilename = sprintf(filename, outputFolder, gate1);
    load(inputFilename)

    Szn(j)=Sz(end);
    Icn(j)=sum(Ic)/length(Ic1);
    ISzn(j)=sum(Isz)/length(Isz0);

    Sxm(j,:)=Sx;
    Sym(j,:)=Sy;
    Szm(j,:)=Sz;
    Icm(j,:)=Ic;
    ISxm(j,:)=Isx;
    ISym(j,:)=Isy;
    ISzm(j,:)=Isz;
end


Sxm=round(Sxm.*10^3)./10^3;
Sym=round(Sym.*10^3)./10^3;
Szm=round(Szm.*10^3)./10^3;

h10=figure(10);
plot(value,Icn)
title(strcat('Ic',axisvariable))
xlabel('exchange coupling')
ylabel('Ic')
saveas(h10,'Ic(exchange coupling)','fig');

h11=figure(11);
plot(value,ISzn)
title('Isz(exchange coupling)')
xlabel('exchange coupling')
ylabel('Isz')
saveas(h11,'Isz(exchange coupling)','fig');

h12=figure(12);
plot(value,Szn)
title('Sz(exchange coupling)')
xlabel('exchange coupling')
ylabel('Sz')
saveas(h12,'Sz(exchange coupling)','fig');


h13=figure(13);
contourf(t, value, Sxm,100,'Linestyle','none')
title('Sx')
xlabel('T(s)')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h13,'Sxcontour','fig');

h14=figure(14);
contourf(t, value, Sym,100,'Linestyle','none')
title('Sy')
xlabel('T(s)')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h14,'Sycontour','fig');

h15=figure(15);
contourf(t, value, Szm,100,'Linestyle','none')
title('Sz')
xlabel('T(s)')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h15,'Szcontour','fig');

h16=figure(16);
contourf(t, value, Icm,100,'Linestyle','none')
title('Ic')
xlabel('T(s)')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h16,'Iccontour','fig');

h17=figure(17);
contourf(t, value, ISxm,100,'Linestyle','none')
title('Isx')
xlabel('T(s)')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h17,'Isxcontour','fig');


h18=figure(18);
contourf(t, value, ISym,100,'Linestyle','none')
title('Isy')
xlabel('T(s)')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h18,'Isycontour','fig');


h19=figure(19);
contourf(t, value, ISzm,100,'Linestyle','none')
title('Isz')
xlabel('T(s)')
ylabel('exchange coupling (meV)')
colormap morgenstemning
colorbar
saveas(h19,'Iszcontour','fig');
