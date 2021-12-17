function [Et,A,Pcarrier] = SSB_Tx(st,CSPRdB,BW,t)
CSPR=10.^(CSPRdB/10);%dB to linear
Ps=mean(abs(st).^2);
Pcarrier=CSPR*Ps;
A=sqrt(Pcarrier).*exp(-1j*pi*2*BW*t);

Et=A+st;
end

