%First Order Hold
clear all;
Fi=12e3;
ti=[0:1/Fi:255/Fi]; %63,127,255 for different number of samples
xi=10*cos(2*pi*ti*4e3);
k=4;
Fs=k*12e3;
t=[0:1/Fs:255/Fi-1/Fs];
L=length(t);
for i=1:1:length(ti)-1
    x(k*(i-1)+1)=xi(i);
    for j=1:1:k-1
        x(k*(i-1)+1+j)=((k-j)*xi(i)+(j)*xi(i+1))/k;
    end
end
X=fft(x);
subplot(2,2,1);
plot(t,x);
title('Time Domain Analysis');
xlabel('time(secs)');
ylabel('Amplitute');



[B,A]=butter(10,Fi/Fs);
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