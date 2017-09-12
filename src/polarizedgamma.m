%Defining a polarized gamma for each lead

gS(1)=pL*g0(1);
gS(2)=pR*g0(2);
g0tot=g0(1)+g0(2);
gStot=gS(1)+gS(2);
g(1)=g0(1)/2*(1+pL)+g0(2)/2*(1+pR);
g(2)=g0(1)/2*(1-pL)+g0(2)/2*(1-pR);