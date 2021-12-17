function [r1,r2,r3,r4,r5,r6,sym_decision] = Rx_QAM64(I,Q)
M=64;
m=sqrt(M);
sym_decision=I+1j*Q;
demod_i= pamdemod(I,m,0,'gray');
dectobits_i = de2bi(demod_i);
demod_q= pamdemod(Q,m,0,'gray');
dectobits_q = de2bi(demod_q);

r1=dectobits_i(:,1);
r2=dectobits_i(:,2);
r3=dectobits_i(:,3);
r4=dectobits_q(:,1);
r5=dectobits_q(:,2);
r6=dectobits_q(:,3);
end