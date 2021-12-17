function Plot(name,SNRb_dB,BER,BER_th,SER,SER_th,r,Sx)
figure(1)
ber=([name,' BER']);
semilogy(SNRb_dB,BER,'g','linewidth',2),grid on,hold on;
semilogy(SNRb_dB,BER_th,'--xr','linewidth',2);
title(ber);
xlabel(' SNR(dB)');
ylabel('Bit Error Rate');
legend('Numeric BER','Theorytical BER'),hold on;

if ~isnan(SER) | ~isnan(SER_th)
figure(2)
ser=([name,' SER']);
semilogy(SNRb_dB,SER,'g','linewidth',2),grid on,hold on;
semilogy(SNRb_dB,SER_th,'--xr','linewidth',2);
title(ser);
xlabel(' SNR(dB)');
ylabel('Symbol Error Rate');
legend('Numeric SER','Theorytical SER')
end

figure(3)
conste=(['Constellation of ',name]);
plot(r,'.m');grid on;hold on
plot(Sx,'*k');
xlabel('Inphase');
ylabel('Quadrature');
title(conste);
legend('With Noise','Without Noise');
end

