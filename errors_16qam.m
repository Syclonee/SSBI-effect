function [error_bit,error_sym] = errors_16qam(r1,r2,r3,r4,sym_decision,b1,b2,b3,b4,Sx)
    error_b1 = (r1~=b1);
    error_b2 = (r2~=b2);
    error_b3 = (r3~=b3);
    error_b4 = (r4~=b4);
    error_bit = error_b1+error_b2+error_b3+error_b4;
    error_sym=(Sx~=sym_decision);
end

