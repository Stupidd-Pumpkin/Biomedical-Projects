%Zero Order Hold
clear all;
Fi=12*10^3;
ti=[0:1/Fi:255*1/Fi]; %63,127,255 for different number of samples
xi=10*cos(2*pi*ti*1e3);
Fs=120e3;
t=[0:1/Fs:255/Fi];
L=length(t);
for i=1:1:255*Fs/Fi+1
    x(i)=xi(ceil(i/10));
end
X=fft(x);
subplot(2,2,1);
plot(t,x);
title('Time Domain Analysis');
xlabel('time(secs)');
ylabel('Amplitute');



[B,A]=butter(5,Fi/Fs);
y=filter(B,A,x);
Y=fft(y);
subplot(2,2,2);
plot(t,y);
title('Time Domain Analysis');
xlabel('time(secs)');
ylabel('Amplitute');


% F=[0:Fs/L:Fs/2-Fs/L,-Fs/2:Fs/L:-Fs/L];
%plot(F,2*abs(X)/L);
F=-Fs/2:Fs/L:Fs/2-Fs/L;
subplot(2,2,3);
plot(F,2*abs(fftshift(X))/L);
title('Magnitute Spectrum');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

subplot(2,2,4);
plot(F,2*abs(fftshift(Y))/L);
title('Magnitute Spectrum');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

