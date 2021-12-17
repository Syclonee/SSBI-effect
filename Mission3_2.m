% Shlomi Uziel & Naama Bendavid 
% 4th year students in Electrical and Electronics Engineering
% Mission 3.2:CDC implementation in RX DSP - Simulation in optic fiber with AWGN channel (continuous  model) - QPSK,QAM-16,QAM-64
% Demonstration to CDC in the DSP area to reduce the length dependence of the fiber

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
dL=0;
L1=100e3+dL; L2=17e3;
beta=0.1;%rolloff factor

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
% optic fiber
yt1=fiber(st,f,Dsmf,L1);

SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    rt=yt1+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    
    sps_new=2; fs_new=sps_new*SR;
    %dispersion Compensation in dsp area with phase revaluation
    CDC = Dispersion_Compensation(r_srrc,sps,sps_new,fs_new,Ddcf,L2);
    r=downsample(CDC,sps_new);
    phi0=mean(angle(sum(r(1:end).*conj(Sx(1:end)))));
    r=r.*exp(-1j*phi0);%phase revaluation and constellation repair
    
    ri=real(r);
    rq=imag(r);
    bld1=(ri>0);
    error_b1=(bld1~=b1);
    bld2=(rq>0);
    error_b2=(bld2~=b2);
    error=error_b1+error_b2;
    R=bld1+1j*bld2;
    BER(k)=sum(error(100:end-100))/(log2(M)*Ns);
    error_sym=(R~=b);
    SER(k)=sum(error_sym(100:end-100))/Ns;
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

SR=60e9;
fs=sps*SR;
%t=0:1/fs:Ns/fs; N=length(t);
Dsmf=17e-6; Ddcf=-100e-6;
dL=0;
L1=100e3+dL; L2=17e3;
beta=0.1;%rolloff factor

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
% optic fiber
yt1=fiber(st,f,Dsmf,L1);

SNRb_dB = 1:10;
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 2.5./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    rt=yt1+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);

    sps_new=2; fs_new=sps_new*SR;
    %dispersion Compensation in dsp area with phase revaluation
    CDC = Dispersion_Compensation(r_srrc,sps,sps_new,fs_new,Ddcf,L2);    
    r=downsample(CDC,sps_new);
    phi0=mean(angle(sum(r(1:end).*conj(Sx(1:end)))));
    r=r.*exp(-1j*phi0);%phase revaluation and constellation repair
    % decision
    [I,Q] = decistion(r,M,Ns);
    % demodulation
    [r1,r2,r3,r4,sym_decision] = Rx_QAM16(I,Q);
    % error calc
    [error_bit,error_sym] = errors_16qam(r1,r2,r3,r4,sym_decision,b1,b2,b3,b4,Sx);
    BER(k)=sum(error_bit(100:end-100))/(log2(M)*Ns);
    SER(k)=sum(error_sym(100:end-100))/(Ns);
end
%plot
BER_th=4*(m-1)*qfunc(sqrt(3*SNRb*log2(M)/(M-1)))/(m*log2(M));
SER_th = 4*(m-1)*qfunc(sqrt(3*SNRs/(M-1)))/(m);
Plot('16-QAM',SNRb_dB,BER,BER_th,SER,SER_th,r,Sx);

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

SR=60e9;
fs=sps*SR;
%t=0:1/fs:Ns/fs; N=length(t);
Dsmf=17e-6; Ddcf=-100e-6;
dL=0;
L1=100e3+dL; L2=17e3;
beta=0.1;%rolloff factor

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
% optic fiber
yt1=fiber(st,f,Dsmf,L1);

SNRb_dB = 1:18;
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 7./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
BER=[];

for k=1:length(SNRb_dB)
    n = sqrt(NO(k)/2)*Noise;
    rt=yt1+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    
    sps_new=2; fs_new=sps_new*SR;
    %dispersion Compensation in dsp area with phase revaluation
    CDC = Dispersion_Compensation(r_srrc,sps,sps_new,fs_new,Ddcf,L2);
    r=downsample(CDC,sps_new);
    phi0=mean(angle(sum(r(1:end).*conj(Sx(1:end)))));
    r=r.*exp(-1j*phi0);%phase revaluation and constellation repair
    % decision
    [I,Q] = decistion(r,M,Ns);
    % demodulation
    [r1,r2,r3,r4,r5,r6,sym_decision] = Rx_QAM64(I,Q);
    % error calc
    [error_bit,error_sym] = errors_64qam(r1,r2,r3,r4,r5,r6,sym_decision,b1,b2,b3,b4,b5,b6,Sx);
    BER(k)=sum(error_bit(100:end-100))/(log2(M)*Ns);
    SER(k)=sum(error_sym(100:end-100))/(Ns);
end
% plot
BER_th=4*(m-1)*qfunc(sqrt(3*SNRb*log2(M)/(M-1)))/(m*log2(M));
SER_th = 4*(m-1)*qfunc(sqrt(3*SNRs/(M-1)))/(m);
Plot('64-QAM',SNRb_dB,BER,BER_th,SER,SER_th,r,Sx);