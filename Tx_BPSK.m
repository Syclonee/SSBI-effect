function [b1,Sx] = Tx_BPSK(Ns)
b1 = rand(1,Ns)>0.5;
Sx = 2*b1-1;
end

