clear all;
Fs=100e3;
t=(0:1/Fs:4097*1/Fs);
L=length(t);
x=10*cos(2*pi*t*2e3)+6*cos(2*pi*t*3e3)+8*cos(2*pi*t*4e3);

figure(1);
subplot(3,1,1);
plot(t,x);
title('Input Signal');
xlabel('time(secs)');
ylabel('Amplitute');

X=fftshift(fft(x));
F=-Fs/2:Fs/L:Fs/2-Fs/L;  
subplot(3,2,5);
plot(F,2*abs(X)/L);
title('Input Spectrum');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

% B is numerator, A is denominator

[n, Wn]=buttord([2.6/50,3.4/50],[2.2/50,3.8/50],1,24); 
[B,A]=butter(n,Wn); %order, fc/fs
y=filter(B,A,x);


% [n, Wn]=cheb1ord([2.6/50,3.4/50],[2.2/50,3.8/50],1,24); 
% [B,A]=cheby1(n,1,Wn);  %order, ripple of passband ,fc/fs
% % % Low pass filter with pass band parameters below fc
% % [B,A]=cheby2(n,1,Wn); %order, ripple of stop band ,fc/fs
% % % Low pass filter with stop band parameters beyond fc
% y=filter(B,A,x);
% 
% % [B,A]=besself(5,10/50); %order,fc/fs 
% % % A Bessel function based filter has a roll over even lesser
% % y=filter(B,A,x);
% 
% [n, Wn]=ellipord([2.6/50,3.4/50],[2.2/50,3.8/50],1,24); 
% [B,A]=ellip(n,1,40,Wn);
% % %order, ripple of pass band, minimum stop band attenuation, fc/fs
% y=filter(B,A,x);

subplot(3,1,2);
plot(t,y);
title('Output Signal');
xlabel('time(secs)');
ylabel('Amplitute');

Y=fftshift(fft(y));
F=-Fs/2:Fs/L:Fs/2-Fs/L; 
subplot(3,2,6);
plot(F,2*abs(Y)/L);
title('Output Spectrum');
xlabel('Frequency(Hz)');
ylabel('Magnitute');