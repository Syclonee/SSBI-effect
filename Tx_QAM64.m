function [b1,b2,b3,b4,b5,b6,Sx] = Tx_QAM64(Ns)
M=64;
m=sqrt(M);
b1=(rand(1,Ns)>0.5)';
b2=(rand(1,Ns)>0.5)';
b3=(rand(1,Ns)>0.5)';
b4=(rand(1,Ns)>0.5)';
b5=(rand(1,Ns)>0.5)';
b6=(rand(1,Ns)>0.5)';

pam_i=bi2de([b3 b2 b1],'left-msb');
pam_q=bi2de([b6 b5 b4],'left-msb');

Si=real(pammod(pam_i,m,0,'gray')');
Sq=real(pammod(pam_q,m,0,'gray')');
Sx=Si+1j*Sq;
end