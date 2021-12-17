function display_constellation_QPSK(IQ)

%for display of QPSK constellation diagram with color density
%Ix=real part of the symbol vector
%Qx=imaginary part of the symbol vector

Ix=real(IQ);
Qx=imag(IQ);
data=[Ix',Qx'];
count=hist2d(data,-2.5:0.05:2.5,-2.5:0.05:2.5);
imagesc(-2.5:0.05:2.5,-2.5:0.05:2.5,count);
colormap([1,1,1;jet(128)]) 
title('QPSK Constellation')
xlabel('In-Phase')
ylabel('Quadrature')
grid on