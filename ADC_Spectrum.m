function [SNR,SFDR] = ADC_Spectrum(Vin,fs)
l=length(Vin)/2;
VIN=abs(fft(Vin));
VIN=VIN(1:l);
VIN_db=20*log10(VIN);
f=linspace(0,fs/2,l);
VIN(1)=0;
figure
plot(f,VIN_db)
[sig_mag,I]=max(VIN);
VIN(I)=0;%removing the actual signal leaves us with the noise
noise_mag=sqrt(sum(VIN.*VIN));
SNR=20*log10(sig_mag)-20*log10(noise_mag);
SFDR=20*log10(sig_mag)-20*log10(max(VIN));


