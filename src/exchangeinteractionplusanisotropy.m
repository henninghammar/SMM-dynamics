%Exchange interaction

jH=1i.*J.^2.*(G0less.*G0great2-G0great.*G0less2-G1xless.*G1xgreat2+G1xgreat.*G1xless2...
    -G1yless.*G1ygreat2+G1ygreat.*G1yless2-G1zless.*G1zgreat2+G1zgreat.*G1zless2);
jDMx=-J.^2.*(G0less.*G1xgreat2-G0great.*G1xless2-G1xless.*G0great2+G1xgreat.*G0less2);
jDMy=-J.^2.*(G0less.*G1ygreat2-G0great.*G1yless2-G1yless.*G0great2+G1ygreat.*G0less2);
jDMz=-J.^2.*(G0less.*G1zgreat2-G0great.*G1zless2-G1zless.*G0great2+G1zgreat.*G0less2);
jIxx=2*1i.*J.^2.*(G1xless.*G1xgreat2-G1xgreat.*G1xless2);
jIyy=2*1i.*J.^2.*(G1yless.*G1ygreat2-G1ygreat.*G1yless2);
jIzz=2*1i.*J.^2.*(G1zless.*G1zgreat2-G1zgreat.*G1zless2)+Dzz;
jIxy=1i.*J.^2.*(G1xless.*G1ygreat2-G1xgreat.*G1yless2+G1yless.*G1xgreat2-G1ygreat.*G1xless2);
jIyx=jIxy;
jIxz=1i.*J.^2.*(G1xless.*G1zgreat2-G1xgreat.*G1zless2+G1zless.*G1xgreat2-G1zgreat.*G1xless2);
jIzx=jIxz;
jIyz=1i.*J.^2.*(G1yless.*G1zgreat2-G1ygreat.*G1zless2+G1zless.*G1ygreat2-G1zgreat.*G1yless2);
jIzy=jIyz;
