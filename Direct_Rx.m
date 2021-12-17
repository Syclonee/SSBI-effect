function [Snt,SNSR,ent] = Direct_Rx(Et,st,Pcarrier,N_iterations,BW,t)
%direct reciver
It=abs(Et).^2;
It=It-mean(It);% dc block
[Snt,SNSR,ent] = iterative_SSBI(It,st,Pcarrier,N_iterations);%filtering
Snt=Snt.*exp(-1j*pi*2*BW*t);%freq shift
end

