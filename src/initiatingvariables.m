%Defining some variables
A=zeros(2,length(w));
B=zeros(2,length(w));
C=zeros(2,2,length(w));
D=zeros(2,2,length(w));
A0=zeros(2,length(w));
A1=zeros(2,length(w));
B0=zeros(2,length(w));
B1=zeros(2,length(w));
C00=zeros(2,length(w));
C01=zeros(2,length(w));
C10=zeros(2,length(w));
C11=zeros(2,length(w));
D00=zeros(2,length(w));
D01=zeros(2,length(w));
D10=zeros(2,length(w));
D11=zeros(2,length(w));
E00=zeros(2,length(w));
E01=zeros(2,length(w));
E10=zeros(2,length(w));
E11=zeros(2,length(w));
F00=zeros(2,length(w));
F01=zeros(2,length(w));
F10=zeros(2,length(w));
F11=zeros(2,length(w));
g0less0=zeros(2,length(w));
g0great0=zeros(2,length(w));
G0less0=zeros(2,length(w));
G0great0=zeros(2,length(w));
g1less0=zeros(2,length(w));
g1great0=zeros(2,length(w));
G1xless0=zeros(2,length(w));
G1xgreat0=zeros(2,length(w));
G1yless0=zeros(2,length(w));
G1ygreat0=zeros(2,length(w));
G1zless0=zeros(2,length(w));
G1zgreat0=zeros(2,length(w));
G0less=zeros(1,tback+1);
G0great=zeros(1,tback+1);
G1xless=zeros(1,tback+1);
G1xgreat=zeros(1,tback+1);
G1yless=zeros(1,tback+1);
G1ygreat=zeros(1,tback+1);
G1zless=zeros(1,tback+1);
G1zgreat=zeros(1,tback+1);
A2=zeros(2,length(w));
B2=zeros(2,length(w));
C2=zeros(2,2,length(w));
D2=zeros(2,2,length(w));
A02=zeros(2,length(w));
A12=zeros(2,length(w));
B02=zeros(2,length(w));
B12=zeros(2,length(w));
C002=zeros(2,length(w));
C012=zeros(2,length(w));
C102=zeros(2,length(w));
C112=zeros(2,length(w));
D002=zeros(2,length(w));
D012=zeros(2,length(w));
D102=zeros(2,length(w));
D112=zeros(2,length(w));
E002=zeros(2,length(w));
E012=zeros(2,length(w));
E102=zeros(2,length(w));
E112=zeros(2,length(w));
F002=zeros(2,length(w));
F012=zeros(2,length(w));
F102=zeros(2,length(w));
F112=zeros(2,length(w));
g0less02=zeros(2,length(w));
g0great02=zeros(2,length(w));
G0less02=zeros(2,length(w));
G0great02=zeros(2,length(w));
g1less02=zeros(2,length(w));
g1great02=zeros(2,length(w));
G1xless02=zeros(2,length(w));
G1xgreat02=zeros(2,length(w));
G1yless02=zeros(2,length(w));
G1ygreat02=zeros(2,length(w));
G1zless02=zeros(2,length(w));
G1zgreat02=zeros(2,length(w));
G0less2=zeros(1,tback+1);
G0great2=zeros(1,tback+1);
G1xless2=zeros(1,tback+1);
G1xgreat2=zeros(1,tback+1);
G1yless2=zeros(1,tback+1);
G1ygreat2=zeros(1,tback+1);
G1zless2=zeros(1,tback+1);
G1zgreat2=zeros(1,tback+1);
K0less=zeros(2,length(w));
K0great=zeros(2,length(w));
Kless=zeros(1,tback+1);
Kgreat=zeros(1,tback+1);
Ic0=zeros(1,length(t));
Ic1=zeros(1,length(t));
Isx0=zeros(1,length(t));
Isx1=zeros(1,length(t));
Isy0=zeros(1,length(t));
Isy1=zeros(1,length(t));
Isz0=zeros(1,length(t));
Isz1=zeros(1,length(t));
dSxvect1=zeros(1,length(t));
dSyvect1=zeros(1,length(t));
dSzvect1=zeros(1,length(t));
dSxvect2=zeros(1,length(t));
dSyvect2=zeros(1,length(t));
dSzvect2=zeros(1,length(t));
dSxvectTot=zeros(1,length(t));
dSyvectTot=zeros(1,length(t));
dSzvectTot=zeros(1,length(t));