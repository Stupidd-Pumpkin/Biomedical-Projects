%Bandpass Sampling/ Under Sampling/ IF Sampling/ Sub-Nyquist Sampling/ 
clear all;
Fc=990e3;
B=20e3;

F1=Fc-B/2;
F2=Fc+B/2;

k=floor(F2/B);
Bp=F2/k;
Fs=2*Bp;

t=[0:1/Fs:255*1/Fs]; %63,127,255 for different number of samples
L=length(t);
x=20*cos(2*pi*990e3*t)+25*cos(2*pi*995e3*t)+40*cos(2*pi*999e3*t)+12*cos(2*pi*985e3*t);
subplot(2,1,1);
plot(t,x);
title('Time Domain Analysis');
xlabel('time(secs)');
ylabel('Amplitute');

X=fft(x);
% F=[0:Fs/L:Fs/2-Fs/L,-Fs/2:Fs/L:-Fs/L];
%plot(F,2*abs(X)/L);
F=-Fs/2:Fs/L:Fs/2-Fs/L;
if (mod(k,2)==1) 
    X=fftshift(X);    
end
subplot(2,2,3);
plot(F,2*abs(X)/L);
title('Magnitute Spectrum');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

subplot(2,2,4);
plot(F,180*angle(X)/pi);
title('Phase Spectrum');
xlabel('Frequency(Hz)');
ylabel('Phase');