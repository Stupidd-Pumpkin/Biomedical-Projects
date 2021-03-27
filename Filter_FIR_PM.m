clear all;
Fs=100e3;
t=(0:1/Fs:4095*1/Fs);
L=length(t);
x=10*cos(2*pi*t*2e3)+6*cos(2*pi*t*3e3)+8*cos(2*pi*t*4e3);
X=fftshift(fft(x));
F=-Fs/2:Fs/L:Fs/2-Fs/L;

[N,Fo,Ao,W]=firpmord([2.2,2.7,3.3,3.8],[0,1,0],[0.1,0.01,0.1],100);
h=firpm(N,Fo,Ao,W);
H=fftshift(fft(h));

y=filtfilt(h,1,x);
Y=fftshift(fft(y));

% figure(1)
% subplot(2,1,1);
% plot(t,x);
% title('Input Signal');
% xlabel('time(secs)');
% ylabel('Amplitute');
% subplot(2,1,2);
% plot(F,2*abs(X)/L);
% title('Magnitute Spectrum');
% xlabel('Frequency(Hz)');
% ylabel('Magnitute');

figure(2)
subplot(2,1,1);
plot(h);
title('Filter');
xlabel('time(secs)');
ylabel('Amplitute');
subplot(2,1,2);
plot(abs(H));
xlabel('Frequency(Hz)');
ylabel('Magnitute');

figure(3)
freqz(h);

% figure(4)
% subplot(2,1,1);
% plot(t,y);
% title('Output Signal');
% xlabel('time(secs)');
% ylabel('Amplitute');
% subplot(2,1,2);
% plot(F,2*abs(Y)/L);
% title('Magnitute Spectrum');
% xlabel('Frequency(Hz)');
% ylabel('Magnitute');
