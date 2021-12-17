% Shlomi Uziel & Naama Bendavid 
% 4th year students in Electrical and Electronics Engineering
% Mission 4:add carrier to transmit our signal - use direct reciver  - QPSK,QAM-16,QAM-64
% Simulation of signal transmission with direct receiver instead of
% coherent receiver in optic system.
% Formation of SSBI effect.

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
t=0:1/fs:(N-1)/fs;%%%%%%%%%%%%%
BW=(1+beta)*SR/2;
CSPRdB=5:1:20;

SNRb_dB = 10;%snr to get ber=10e-3 
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
n = sqrt(NO/2)*Noise;

for k=1:length(CSPRdB)
    [Et,A,Pcarrier] = SSB_Tx(st,CSPRdB(k),BW,t);%SSB Tx
    Et=fiber(Et,f,Dsmf,L1);%fiber
    rt=Et+n;

    %Direct reciver
    N_iterations=1;
    It=abs(rt).^2;
    It=It-mean(It);% dc block
    [Snt] = iterative_SSBI(It,st,Pcarrier,N_iterations);%filtering
    Snt=Snt.*exp(-1j*pi*2*BW*t);%freq shift

    r_srrc=conv(Snt,h);
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
    BER(k)=sum(error)/(log2(M)*Ns);
end
% plot
figure(1)
semilogy(CSPRdB,BER,'g','linewidth',2),grid on,hold on;
title('QPSK BER');
xlabel('CSPR(dB)');ylabel('Bit Error Rate');
figure(2)
plot(r,'.m');grid on;hold on
plot(Sx,'*k');
xlabel('Inphase');ylabel('Quadrature');
title('Constellation of QPSK');

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

s1=upsample(Sx,sps);
h=rcosdesign(beta,20,sps,'sqrt');
st=conv(s1,h);
st=st(sps*10+1:end-sps*10);

N=length(st);
f=fs*(-0.5:1/N:0.5-1/N);
t=0:1/fs:(N-1)/fs;%%%%%%%%%%%%%
BW=(1+beta)*SR/2;
CSPRdB=10:1:20;

SNRb_dB = 15;%snr to get ber=10e-3 
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
    N_iterations=1;
    It=abs(rt).^2;
    It=It-mean(It);% dc block
    [Snt] = iterative_SSBI(It,st,Pcarrier,N_iterations);%filtering
    Snt=Snt.*exp(-1j*pi*2*BW*t);%freq shift
    
    r_srrc=conv(Snt,h);
    r_srrc=r_srrc(sps*10+1:end-sps*10);
    %dsp
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
    BER(k)=sum(error_bit)/(log2(M)*Ns);
    %SER(k)=sum(error_sym(100:end-100))/(Ns);
end
% plot
figure(1)
semilogy(CSPRdB,BER,'g','linewidth',2),grid on,hold on;
title('16-QAM BER');
xlabel('CSPR(dB)');ylabel('Bit Error Rate');
figure(2)
plot(r,'.m');grid on;hold on
plot(Sx,'*k');
xlabel('Inphase');ylabel('Quadrature');
title('Constellation of 16-QAM');

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
t=0:1/fs:(N-1)/fs;%%%%%%%%%%%%%
BW=(1+beta)*SR/2;
CSPRdB=10:1:20;

SNRb_dB = 25;%snr to get ber=10e-3 
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
    N_iterations=1;
    It=abs(rt).^2;
    It=It-mean(It);% dc block
    [Snt] = iterative_SSBI(It,st,Pcarrier,N_iterations);%filtering
    Snt=Snt.*exp(-1j*pi*2*BW*t);%freq shift

    
    r_srrc=conv(Snt,h);
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
    BER(k)=sum(error_bit)/(log2(M)*Ns);
    %SER(k)=sum(error_sym(100:end-100))/(Ns);
end
% plot
figure(1)
semilogy(CSPRdB,BER,'g','linewidth',2),grid on,hold on;
title('64-QAM BER');
xlabel('CSPR(dB)');ylabel('Bit Error Rate');
figure(2)
plot(r,'.m');grid on;hold on
plot(Sx,'*k');
xlabel('Inphase');ylabel('Quadrature');
title('Constellation of 64-QAM');

%% QPSK CSPR constellation
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
t=0:1/fs:(N-1)/fs;%%%%%%%%%%%%%
BW=(1+beta)*SR/2;
CSPRdB=5:5:20;

SNRb_dB = 12;%snr to get ber=10e-3 
SNRb = 10.^(SNRb_dB/10);
NO = 1./ SNRb;
Noise=randn(1,length(st))+1j*(randn(1,length(st)));
n = sqrt(NO/2)*Noise;

for k=1:length(CSPRdB)
    [Et,A,Pcarrier] = SSB_Tx(st,CSPRdB(k),BW,t);%SSB Tx
    Et=fiber(Et,f,Dsmf,L1);%fiber
    rt=Et+n;

    %Direct reciver
    N_iterations=1;
    It=abs(rt).^2;
    It=It-mean(It);% dc block
    [Snt] = iterative_SSBI(It,st,Pcarrier,N_iterations);%filtering
    Snt=Snt.*exp(-1j*pi*2*BW*t);%freq shift


    r_srrc=conv(Snt,h);
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
    BER(k)=sum(error)/(log2(M)*Ns);
    figure(2)
    plot(r,'x');grid on;hold on
    leg=append('CSPR=',string(CSPRdB));
    legend(leg);hold on
end
% plot
figure(2)
plot(Sx,'*k');
xlabel('Inphase');ylabel('Quadrature');
title('Constellation of QPSK');

