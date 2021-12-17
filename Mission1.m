% Shlomi Uziel & Naama Bendavid 
% 4th year students in Electrical and Electronics Engineering
% Mission 1: Simulation in AWGN channel (continuous  model) - BPSK,QPSK,QAM-16,QAM-64
% compre simulation to theorytical BER & SER

%% BPSK
clc
clear
close all

Ns = 1e6;
[b1,Sx] = Tx_BPSK(Ns);
sps=5;

s1=upsample(Sx,5);
h=rcosdesign(0.1,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);
SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st));
BER=[];

for k=1:length(SNRb_dB)
    n=sqrt(NO(k)/2)*Noise;
    rt=st+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    r=downsample(r_srrc,sps);
    
    bld=(r>0);
    error=(bld~=b1);
    BER(k)=sum(error)/Ns;
end
BER_th = qfunc(sqrt(2*SNRb));
Plot('BPSK',SNRb_dB,BER,BER_th,NaN ,NaN ,r,Sx)

%% QPSK
clc
clear
close all

Ns = 5e5;
M=4;
m=sqrt(M);
[b1,b2,b,Sx] = Tx_QPSK(Ns);
sps=5;

s1=upsample(Sx,5);
h=rcosdesign(0.1,20,sps,'sqrt');
st=conv(s1,h);
st=st(51:end-50);
SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    rt=st+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    r=downsample(r_srrc,sps);
    
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

%% 16-QAM
clear
clc
close all

Ns=1e5;
M=16;
m=sqrt(M);
% rand bits and make them dec for 4pam mod
[b1,b2,b3,b4,Sx] = Tx_QAM16(Ns);
sps=5;

s1=upsample(Sx,5);
h=rcosdesign(0.1,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 2.5./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    rt=st+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    r=downsample(r_srrc,sps);
    
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


figure('name' ,'16-QAM')
plot(Sx,'og','linewidth',10);grid on;
xlabel('Inphase');ylabel('Quadrature');
title('16-QAM');
xlim([-4,4]);ylim([-4,4])

%% 64-QAM
clear
clc
close all

Ns=1e5;
M=64;
m=sqrt(M);
% rand bits and make them dec for 6pam mod
[b1,b2,b3,b4,b5,b6,Sx] = Tx_QAM64(Ns);
sps=5;

s1=upsample(Sx,5);
h=rcosdesign(0.1,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

SNRb_dB = 1:18;
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 7./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    rt=st+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    r=downsample(r_srrc,sps);
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