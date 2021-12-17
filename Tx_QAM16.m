function [b1,b2,b3,b4,Sx] = Tx_QAM16(Ns)
M=16;
m=sqrt(M);
b1=(rand(1,Ns)>0.5)';
b2=(rand(1,Ns)>0.5)';
b3=(rand(1,Ns)>0.5)';
b4=(rand(1,Ns)>0.5)';

pam_i=bi2de([b2 b1],'left-msb');
pam_q=bi2de([b4 b3],'left-msb');

Si=real(pammod(pam_i,m,0,'gray')');
Sq=real(pammod(pam_q,m,0,'gray')');
Sx=Si+1j*Sq;
end