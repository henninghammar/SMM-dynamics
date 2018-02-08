
figure(1)
plot(eVvector, Icbias)
ylabel('Charge current')
xlabel('Bias voltage (1/Gamma)')

figure(2)
hold on
plot(eVvector, Isxbias)
plot(eVvector, Isybias)
plot(eVvector, Iszbias)
legend('Isx','Isy','Isz')
ylabel('Spin current')
xlabel('Bias voltage (1/Gamma)')

figure(3)
plot(eVvector, Iebias)
ylabel('Energy current')
xlabel('Bias voltage (1/Gamma)')

figure(4)
plot(eVvector, Iqbias)
ylabel('Heat current')
xlabel('Bias voltage (1/Gamma)')
%
% figure(5)
% contourf(w, eVvector, DOSbias, 100,'Linestyle','none')
% title('Density of states')
% ylabel('Bias voltage (1/Gamma)')
% colorbar
%
% figure(6)
% contourf(w, eVvector, MDOSzbias, 100,'Linestyle','none')
% title('Spin density of states (z)')
% ylabel('Bias voltage (1/Gamma)')
% colorbar

figure(10)
plot(epsvector, Icgate)
ylabel('Charge current')
xlabel('Gate voltage (1/Gamma)')

figure(11)
hold on
plot(epsvector, Isxgate)
plot(epsvector, Isygate)
plot(epsvector, Iszgate)
legend('Isx','Isy','Isz')
ylabel('Spin current')
xlabel('Gate voltage (1/Gamma)')

figure(12)
plot(epsvector, Iegate)
ylabel('Energy current')
xlabel('Gate voltage (1/Gamma)')

figure(13)
plot(epsvector, Iqgate)
ylabel('Heat current')
xlabel('Gate voltage (1/Gamma)')
%
% figure(14)
% contourf(w, epsvector, DOSgate, 100,'Linestyle','none')
% title('Density of states')
% ylabel('Gate voltage (1/Gamma)')
% colorbar
%
% figure(15)
% contourf(w, epsvector, MDOSzgate, 100,'Linestyle','none')
% title('Spin density of states (z)')
% ylabel('Gate voltage (1/Gamma)')
% colorbar

figure(20)
plot(Jvector, IcJ)
ylabel('Charge current')
xlabel('Exchange coupling (1/Gamma)')

figure(21)
hold on
plot(Jvector, IsxJ)
plot(Jvector, IsyJ)
plot(Jvector, IszJ)
legend('Isx','Isy','Isz')
ylabel('Spin current')
xlabel('Exchange coupling (1/Gamma)')

figure(22)
plot(Jvector, IeJ)
ylabel('Energy current')
xlabel('Exchange coupling (1/Gamma)')

figure(23)
plot(Jvector, IqJ)
ylabel('Heat current')
xlabel('Exchange coupling (1/Gamma)')
%
% figure(24)
% contourf(w, Jvector, DOSJ, 100,'Linestyle','none')
% title('Density of states')
% ylabel('Exchange coupling (1/Gamma)')
% colorbar
%
% figure(25)
% contourf(w, Jvector, MDOSzJ, 100,'Linestyle','none')
% title('Spin density of states (z)')
% ylabel('Exchange coupling (1/Gamma)')
% colorbar

figure(30)
plot(tempvector, real(Ictemp))
ylabel('Charge current')
xlabel('Temperature T_R (K) (T_L = 1 K)')

figure(31)
hold on
plot(tempvector, real(Isxtemp))
plot(tempvector, real(Isytemp))
plot(tempvector, real(Isztemp))
legend('Isx','Isy','Isz')
ylabel('Spin current')
xlabel('Temperature T_R (K) (T_L = 1 K)')

figure(32)
plot(tempvector, real(Ietemp))
ylabel('Energy current')
xlabel('Temperature T_R (K) (T_L = 1 K)')

figure(33)
plot(tempvector, real(Iqtemp))
ylabel('Heat current')
xlabel('Temperature T_R (K) (T_L = 1 K)')
%
% figure(34)
% contourf(w, tempvector, real(DOStemp), 100,'Linestyle','none')
% title('Density of states')
% ylabel('Temperature T_R (K) (T_L = 1 K)')
% colorbar
%
% figure(35)
% contourf(w, tempvector, real(MDOSztemp), 100,'Linestyle','none')
% title('Spin density of states (z)')
% ylabel('Temperature T_R (K) (T_L = 1 K)')
% colorbar
