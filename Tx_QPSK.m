function [b1,b2,b,Sx] = Tx_QPSK(Ns)
b1 = rand(1,Ns)>0.5;
b2 = rand(1,Ns)>0.5;
b = b1+1j*b2;

Si = 2*b1-1;
Sq = 2*b2-1;
Sx = Si+1j*Sq;
end

