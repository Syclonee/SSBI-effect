function [error_bit,error_sym] = errors_64qam(r1,r2,r3,r4,r5,r6,sym_decision,b1,b2,b3,b4,b5,b6,Sx)
    error_b1 = (r1~=b1);
    error_b2 = (r2~=b2);
    error_b3 = (r3~=b3);
    error_b4 = (r4~=b4);
    error_b5 = (r5~=b5);
    error_b6 = (r6~=b6);
    error_bit = error_b1+error_b2+error_b3+error_b4+error_b5+error_b6;
    error_sym=(Sx~=sym_decision);
end

