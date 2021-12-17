%% QPSK
clc
clear
close all

Ns = 5e5;
M=4;
m=sqrt(M);
[b1,b2,b,Sx] = Tx_QPSK(Ns);
sps=5;

SR=60e9;
fs=sps*SR;
%t=0:1/fs:Ns/fs; N=length(t);
Dsmf=17e-6; Ddcf=-100e-6;
L1=100e3; L2=17e3;
beta=0.1;%rolloff factor

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
% optic fiber
yt1=fiber(st,f,Dsmf,L1);
yt2=fiber(yt1,f,Ddcf,L2);

SNRb_dB = 1:20;
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    rt=yt2+n;
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


plot(r,'.r');grid on;hold on
plot(Sx,'xb','linewidth',4);
xlabel('Inphase');ylabel('Quadrature');
legend('Symbol with Noise','Symbol without Noise')
title(['SNRb = ',int2str(SNRb_dB(end)),'dB']);
xlim([-2,2]);ylim([-2,2]);