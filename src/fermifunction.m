%Defining a fermi-function for each lead

fermi=zeros(2,length(w));
beta(1)=1/(kB*T(1));
beta(2)=1/(kB*T(2));
fermi(1,:)=1./(1+exp(beta(1).*w));
fermi(2,:)=1./(1+exp(beta(2).*w));