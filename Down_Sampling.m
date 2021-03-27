%Down Sampling
clear all;
Fi=50*10^3;
ti=[0:1/Fi:255/Fi]; %63,127,255 for different number of samples
xi=10*cos(2*pi*ti*10^3)+6*cos(2*pi*ti*2*10^3)+2*cos(2*pi*ti*4*10^3);
k=5;
Fs=Fi/k;
t=[0:1/Fs:floor(length(ti)/k)/Fs];
L=length(t);
for i=1:1:length(t)
    x(i)=xi(k*i-(k-1));
end
subplot(2,1,1);
plot(t,x);
title('Time Domain Analysis');
xlabel('time(secs)');
ylabel('Amplitute');

Xi=fft(xi);
X=fft(x);
% F=[0:Fs/L:Fs/2-Fs/L,-Fs/2:Fs/L:-Fs/L];
%plot(F,2*abs(X)/L);
F=-Fs/2:Fs/L:Fs/2-Fs/L;
Y=fftshift(X);    
subplot(2,2,3);
plot(F,2*abs(Y)/L);
title('Magnitute Spectrum');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

subplot(2,2,4);
plot(F,180*angle(Y)/pi);
title('Phase Spectrum');
xlabel('Frequency(Hz)');
ylabel('Phase');