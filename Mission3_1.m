% Shlomi Uziel & Naama Bendavid 
% 4th year students in Electrical and Electronics Engineering
% Mission 3.1:DCF Limitation - Simulation in optic fiber with AWGN channel (continuous  model) -QPSK,QAM-16,QAM-64
% Examine how the BER varies in constant SNRB and variable fiber length.

%% QPSK
clc
clear
close all

Ns = 1e5;
M=4;
m=sqrt(M);
[b1,b2,b,Sx] = Tx_QPSK(Ns);
sps=5;

figu=['g','r','c'];
SR=[15e9,30e9,60e9];

for x=1:length(SR)
fs=sps*SR(x);
%t=0:1/fs:Ns/fs; N=length(t);
Dsmf=17e-6; Ddcf=-100e-6;
dL=-1e3:100:1e3;
L1=100e3; L2=17e3;
beta=0.1;%rolloff factor

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);

SNRb_dB = 7;%ber~1e-3
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));%
n = sqrt(NO/2)*Noise;
BER=[];

for k=1:length(dL)
    L1=100e3;
    L1=L1+dL(k);
    % optic fiber
    yt1=fiber(st,f,Dsmf,L1);
    yt2=fiber(yt1,f,Ddcf,L2);
    
    rt=yt2+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    r=downsample(r_srrc,sps);
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
    BER(k)=sum(error)/(log2(M)*Ns);
    %error_sym=(R~=b);
    %SER(k)=sum(error_sym)/Ns;
end
%plot
semilogy(dL,BER,figu(x),'linewidth',1.5),grid on,hold on;
name=append('SR=',string(SR*1e-9),'GBaund');
title('QPSK BER');legend(name),hold on;
xlabel('\DeltaL(m)'); ylabel('Bit Error Rate');
end

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

figu=['g','r','c'];
SR=[15e9,30e9,60e9];

for x=1:length(SR)
fs=sps*SR(x);
%t=0:1/fs:Ns/fs; N=length(t);
Dsmf=17e-6; Ddcf=-100e-6;
dL=-1e3:100:1e3;
L1=100e3; L2=17e3;
beta=0.1;%rolloff factor

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);

SNRb_dB = 11;%ber~1e-3
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 2.5./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
n = sqrt(NO/2)*Noise;
BER=[];

for k=1:length(dL)
    L1=100e3;
    L1=L1+dL(k);
    % optic fiber
    yt1=fiber(st,f,Dsmf,L1);
    yt2=fiber(yt1,f,Ddcf,L2);

    rt=yt2+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    r=downsample(r_srrc,sps);
    phi0=mean(angle(sum(r(1:end).*conj(Sx(1:end)))));
    r=r.*exp(-1j*phi0);%phase revaluation and constellation repair
    
    % decision
    [I,Q] = decistion(r,M,Ns);
    % demodulation
    [r1,r2,r3,r4,sym_decision] = Rx_QAM16(I,Q);
    % error calc
    [error_bit,error_sym] = errors_16qam(r1,r2,r3,r4,sym_decision,b1,b2,b3,b4,Sx);
    BER(k)=sum(error_bit)/(log2(M)*Ns);
    %SER(k)=sum(error_sym)/(Ns);
end
%plot
semilogy(dL,BER,figu(x),'linewidth',1.5),grid on,hold on;
name=append('SR=',string(SR*1e-9),'GBaund');
title('16-QAM BER');legend(name),hold on;
xlabel('\DeltaL(m)'); ylabel('Bit Error Rate');
end

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

figu=['g','r','c'];
SR=[15e9,30e9,60e9];
for x=1:length(SR)
fs=sps*SR(x);
%t=0:1/fs:Ns/fs; N=length(t);
Dsmf=17e-6; Ddcf=-100e-6;
dL=-1e3:100:1e3;
L1=100e3; L2=17e3;
beta=0.1;%rolloff factor

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);

SNRb_dB = 15;%ber~1e-3
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 7./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
n = sqrt(NO/2)*Noise;
BER=[];

for k=1:length(dL)
    L1=100e3;
    L1=L1+dL(k);
    % optic fiber
    yt1=fiber(st,f,Dsmf,L1);
    yt2=fiber(yt1,f,Ddcf,L2);
    
    rt=yt2+n;
    r_srrc=conv(rt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    r=downsample(r_srrc,sps);
    phi0=mean(angle(sum(r(1:end).*conj(Sx(1:end)))));
    r=r.*exp(-1j*phi0);%phase revaluation and constellation repair
    % decision
    [I,Q] = decistion(r,M,Ns);

    % demodulation
    [r1,r2,r3,r4,r5,r6,sym_decision] = Rx_QAM64(I,Q);
    % error calc
    [error_bit,error_sym] = errors_64qam(r1,r2,r3,r4,r5,r6,sym_decision,b1,b2,b3,b4,b5,b6,Sx);
    BER(k)=sum(error_bit)/(log2(M)*Ns);
    %SER(k)=sum(error_sym)/(Ns);
end
%plot
semilogy(dL,BER,figu(x),'linewidth',1.5),grid on,hold on;
name=append('SR=',string(SR*1e-9),'GBaund');
title('64-QAM BER');legend(name),hold on;
xlabel('\DeltaL(m)'); ylabel('Bit Error Rate');
end
