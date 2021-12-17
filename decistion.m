function [I,Q] = decistion(r,M,Ns)
    m=sqrt(M);
    low_sym=-m+1;
    ri=real(r);
    rq=imag(r);
    I=zeros(1,Ns);
    Q=zeros(1,Ns);
    % decision
    %low_sym=-3;
    I(ri<=low_sym+1)=low_sym;
    Q(rq<=low_sym+1)=low_sym;
    for x=1:m-1
      I(ri>=low_sym+1)=low_sym+2;
       Q(rq>=low_sym+1)=low_sym+2;
       low_sym=low_sym+2; 
    end
end

