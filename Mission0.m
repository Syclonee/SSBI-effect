% Shlomi Uziel & Naama Bendavid 
% 4th year students in Electrical and Electronics Engineering
% Mission 0: Simulation in AWGN channel (discrete  model) - BPSK,QPSK,QAM-16,QAM-64
% compre simulation to theorytical BER & SER

%% BPSK
clear
clc
close all

Ns = 1e6;
[b1,Sx] = Tx_BPSK(Ns);
SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,Ns);
BER=[];

for k=1:length(SNRb_dB)
    n=sqrt(NO(k)/2)*Noise;
    r=Sx+n;
    bld=(r>0);
    error=(bld~=b1);
    BER(k)=sum(error)/Ns;
end
% plot
BER_th = qfunc(sqrt(2*SNRb));
Plot('BPSK',SNRb_dB,BER,BER_th,NaN ,NaN ,r,Sx)

%% QPSK
clear
clc
close all

Ns = 1e6;
M=4;
[b1,b2,b,Sx] = Tx_QPSK(Ns);
SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=(randn(1,Ns)+1j*(randn(1,Ns)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    r=Sx+n;
    ri=real(r);
    rq=imag(r);
    bld1=(ri>0);
    error_b1=(bld1~=b1);
    bld2=(rq>0);
    error_b2=(bld2~=b2);
    error=error_b1+error_b2;
    R=bld1+1j*bld2;
    BER(k)=sum(error)/(log2(M)*Ns);
    error_sym=(R~=b);
    SER(k)=sum(error_sym)/Ns;
end
% plot
BER_th = qfunc(sqrt(2*SNRb));
SER_th = 2.*qfunc(sqrt(2*SNRb))-(qfunc(sqrt(2*SNRb))).^2;
Plot('QPSK',SNRb_dB,BER,BER_th,SER,SER_th,r,Sx);

%% QAM-16
clear
clc
close all

Ns=1e5;
M=16;
m=sqrt(M);
% rand bits and make them dec for 4pam mod
[b1,b2,b3,b4,Sx] = Tx_QAM16(Ns);

SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 2.5./ SNRb;
Noise=(randn(1,Ns)+1j*(randn(1,Ns)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    r = Sx+n;
    % decision
    [I,Q] = decistion(r,M,Ns);
    % demodulation
    [r1,r2,r3,r4,sym_decision] = Rx_QAM16(I,Q);
    % error calc
    [error_bit,error_sym] = errors_16qam(r1,r2,r3,r4,sym_decision,b1,b2,b3,b4,Sx);
    BER(k)=sum(error_bit)/(log2(M)*Ns);
    SER(k)=sum(error_sym)/(Ns);
end
%plot
BER_th=4*(m-1)*qfunc(sqrt(3*SNRb*log2(M)/(M-1)))/(m*log2(M));
SER_th = 4*(m-1)*qfunc(sqrt(3*SNRs/(M-1)))/(m);
Plot('16-QAM',SNRb_dB,BER,BER_th,SER,SER_th,r,Sx);

%% QAM-64
clear
clc
close all

Ns=1e5;
M=64;
m=sqrt(M);
% rand bits and make them dec for 6pam mod
[b1,b2,b3,b4,b5,b6,Sx] = Tx_QAM64(Ns);

SNRb_dB = 1:18;
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 7./ SNRb;
Noise=(randn(1,Ns)+1j*(randn(1,Ns)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    r = Sx+n;

    % decision
    [I,Q] = decistion(r,M,Ns);

    % demodulation
    [r1,r2,r3,r4,r5,r6,sym_decision] = Rx_QAM64(I,Q);
    % error calc
    [error_bit,error_sym] = errors_64qam(r1,r2,r3,r4,r5,r6,sym_decision,b1,b2,b3,b4,b5,b6,Sx);
    BER(k)=sum(error_bit)/(log2(M)*Ns);
    SER(k)=sum(error_sym)/(Ns);
end
% plot
BER_th=4*(m-1)*qfunc(sqrt(3*SNRb*log2(M)/(M-1)))/(m*log2(M));
SER_th = 4*(m-1)*qfunc(sqrt(3*SNRs/(M-1)))/(m);
Plot('64-QAM',SNRb_dB,BER,BER_th,SER,SER_th,r,Sx);