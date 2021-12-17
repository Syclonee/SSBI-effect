function [Snt] = iterative_SSBI(It,st,Pcarrier,N_iterations)
S0t=0.5*hilbert(It)./sqrt(Pcarrier);
Dn=0;

for n=1:N_iterations
Snt=S0t-Dn;
Dn=0.5*hilbert(abs(Snt).^2)./sqrt(Pcarrier);
%ent=Snt-st;
%SNSR(n)=mean(abs(ent).^2)/mean(abs(st).^2);
end
end

