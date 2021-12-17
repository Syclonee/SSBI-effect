function y = fiber(xt,f,D,L)
lamda=1550e-9;
c=3e8;
beta2=-(lamda^2)*D/(2*pi*c);
N=length(xt);
%fourier transform
Xf=fftshift(fft(xt))*2/N;
Sf=Xf.*exp(0.5*1j*beta2*L*(2*pi*f).^2);
y=ifft(ifftshift(Sf))*N/2;
end

