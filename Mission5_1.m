% Shlomi Uziel & Naama Bendavid 
% 4th year students in Electrical and Electronics Engineering

% Mission 5.1:downsample the signal after the direct Rx and down again after - QPSK,QAM-16,QAM-64

%% QPSK
clc
clear
close all

Ns = 1e5;
M=4;
m=sqrt(M);
[b1,b2,b,Sx] = Tx_QPSK(Ns);
sps=5;

SR=60e9;
fs=sps*SR;
Dsmf=17e-6; Ddcf=-100e-6;
dL=0;
L1=100e3+dL; L2=17e3;
beta=0.1;%rolloff factor

sps_2=2;
s1=upsample(Sx,sps_2);
h=rcosdesign(beta,20,sps_2,'sqrt');
st=conv(s1,h);
st=st(sps_2*10+1:end-sps_2*10);%sps=2

st=resample(st,sps,sps_2);%sps=5

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
t=0:1/fs:(N-1)/fs;
BW=(1+beta)*SR/2;
CSPRdB=5:1:20;

SNRb_dB = 6;%snr to get ber=10e-3 
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
n = sqrt(NO/2)*Noise;

for k=1:length(CSPRdB)
    [Et,A,Pcarrier] = SSB_Tx(st,CSPRdB(k),BW,t);%SSB Tx
    Et=fiber(Et,f,Dsmf,L1);%fiber
    rt=Et+n;

    %Direct reciver
    It=abs(rt).^2;
    It=It-mean(It);%dc block
    
    sps_3=3; fs_3=sps_3*SR;
    It1=resample(It,sps_3,sps);%downsampling to sps=3
    N_=length(It1);
    N_iterations=1;
    [Snt] = iterative_SSBI(It1,st,Pcarrier,N_iterations);%filtering
    t1=0:1/fs_3:(N_-1)/fs_3;
    Snt=Snt.*exp(-1j*pi*2*BW*t1);%freq shift
    
    fs_2=sps_2*SR;
    %dispersion Compensation in dsp area with downsample to sps=2
    CDC = Dispersion_Compensation(Snt,sps_3,sps_2,fs_2,Ddcf,L2);
    
    r_srrc=conv(CDC,h);
    r_srrc=r_srrc(sps_2*10+1:end-sps_2*10);
    
    r=downsample(r_srrc,sps_2);%%%%downsampling to sps=1
    
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
end
% plot
figure(1)
semilogy(CSPRdB,BER,'g','linewidth',2),grid on,hold on;
title('QPSK BER');
xlabel('CSPR(dB)');ylabel('Bit Error Rate');
figure(2)
display_constellation_QPSK(r)

%% 16-QAM
clear
close all
clc 

Ns=1e5;
M=16;
m=sqrt(M);
% rand bits and make them dec for 4pam mod
[b1,b2,b3,b4,Sx] = Tx_QAM16(Ns);
sps=5;

SR=60e9;
fs=sps*SR;
Dsmf=17e-6; Ddcf=-100e-6;
dL=0;
L1=100e3+dL; L2=17e3;
beta=0.1; %rolloff factor

sps_2=2;
s1=upsample(Sx,sps_2);
h=rcosdesign(beta,20,sps_2,'sqrt');
st=conv(s1,h);
st=st(sps_2*10+1:end-sps_2*10);%sps=2

st=resample(st,sps,sps_2);%sps=5

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
t=0:1/fs:(N-1)/fs;
BW=(1+beta)*SR/2;
CSPRdB=5:1:20;

SNRb_dB = 10;%snr to get ber=10e-3 
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 2.5./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
n = sqrt(NO/2)*Noise;

for k=1:length(CSPRdB)
    [Et,A,Pcarrier] = SSB_Tx(st,CSPRdB(k),BW,t);%SSB Tx
    Et=fiber(Et,f,Dsmf,L1);%fiber
    rt=Et+n;

    %Direct reciver
    It=abs(rt).^2;
    It=It-mean(It);%dc block
    
    sps_3=3; fs_3=sps_3*SR;
    It1=resample(It,sps_3,sps);%downsampling to sps=3
    N_=length(It1);
    N_iterations=1;
    [Snt] = iterative_SSBI(It1,st,Pcarrier,N_iterations);%filtering
    t1=0:1/fs_3:(N_-1)/fs_3;
    Snt=Snt.*exp(-1j*pi*2*BW*t1);%freq shift
    
    fs_2=sps_2*SR;
    %dispersion Compensation in dsp area with downsample to sps=2
    CDC = Dispersion_Compensation(Snt,sps_3,sps_2,fs_2,Ddcf,L2);
    
    r_srrc=conv(CDC,h);
    r_srrc=r_srrc(sps_2*10+1:end-sps_2*10);
    
    r=downsample(r_srrc,sps_2);%downsampling to sps=1
    
    phi0=mean(angle(sum(r(1:end).*conj(Sx(1:end)))));
    r=r.*exp(-1j*phi0);%phase revaluation and constellation repair
    
    % decision
    [I,Q] = decistion(r,M,Ns);
    % demodulation
    [r1,r2,r3,r4,sym_decision] = Rx_QAM16(I,Q);
    % error calc
    [error_bit,error_sym] = errors_16qam(r1,r2,r3,r4,sym_decision,b1,b2,b3,b4,Sx);
    BER(k)=sum(error_bit)/(log2(M)*Ns);

end
% plot
figure(1)
semilogy(CSPRdB,BER,'g','linewidth',2),grid on,hold on;
title('16-QAM BER');
xlabel('CSPR(dB)');ylabel('Bit Error Rate');
figure(2)
display_constellation_16QAM(r)

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
Dsmf=17e-6; Ddcf=-100e-6;
dL=0;
L1=100e3+dL; L2=17e3;
beta=0.1;%rolloff factor

sps_2=2;
s1=upsample(Sx,sps_2);
h=rcosdesign(beta,20,sps_2,'sqrt');
st=conv(s1,h);
st=st(sps_2*10+1:end-sps_2*10);%sps=2

st=resample(st,sps,sps_2);%sps=5

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
t=0:1/fs:(N-1)/fs;
BW=(1+beta)*SR/2;
CSPRdB=10:1:20;

SNRb_dB = 20;%snr to get ber=10e-3 
SNRb = 10.^(SNRb_dB/10);
SNRs = log2(M)*SNRb;
NO = 7./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
n = sqrt(NO/2)*Noise;

for k=1:length(CSPRdB)
    [Et,A,Pcarrier] = SSB_Tx(st,CSPRdB(k),BW,t);%SSB Tx
    Et=fiber(Et,f,Dsmf,L1);%fiber
    rt=Et+n;

    %Direct reciver
    It=abs(rt).^2;
    It=It-mean(It);%dc block
    
    sps_3=3; fs_3=sps_3*SR;
    It1=resample(It,sps_3,sps);%downsampling to sps=3
    N_=length(It1);
    N_iterations=1;
    [Snt] = iterative_SSBI(It1,st,Pcarrier,N_iterations);%filtering
    t1=0:1/fs_3:(N_-1)/fs_3;
    Snt=Snt.*exp(-1j*pi*2*BW*t1);%freq shift
    
    fs_2=sps_2*SR;
    %dispersion Compensation in dsp area with downsample to sps=2
    CDC = Dispersion_Compensation(Snt,sps_3,sps_2,fs_2,Ddcf,L2);
    
    r_srrc=conv(CDC,h);
    r_srrc=r_srrc(sps_2*10+1:end-sps_2*10);
    
    r=downsample(r_srrc,sps_2);%downsampling to sps=1
 
    phi0=mean(angle(sum(r(1:end).*conj(Sx(1:end)))));
    r=r.*exp(-1j*phi0);%phase revaluation and constellation repair
    % decision
    [I,Q] = decistion(r,M,Ns);
    % demodulation
    [r1,r2,r3,r4,r5,r6,sym_decision] = Rx_QAM64(I,Q);
    % error calc
    [error_bit,error_sym] = errors_64qam(r1,r2,r3,r4,r5,r6,sym_decision,b1,b2,b3,b4,b5,b6,Sx);
    BER(k)=sum(error_bit)/(log2(M)*Ns);
end
% plot
figure(1)
semilogy(CSPRdB,BER,'g','linewidth',2),grid on,hold on;
title('64-QAM BER');
xlabel('CSPR(dB)');ylabel('Bit Error Rate');
figure(2)
display_constellation_64QAM(r)