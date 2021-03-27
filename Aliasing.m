%Aliasing
Fs=4e3; %4,5,8
t=[0:1/Fs:255*1/Fs]; %63,127,255 for different number of samples
L=length(t);
x=10*cos(2*pi*t*1e3)+6*cos(2*pi*t*2e3)+2*cos(2*pi*t*4e3);
subplot(2,1,1);
plot(t,x);
title('Time Domain Analysis');
xlabel('time(secs)');
ylabel('Amplitute');

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