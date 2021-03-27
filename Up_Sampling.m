%Up Sampling
clear all;
Fi=12*10^3;
ti=[0:1/Fi:255/Fi]; %63,127,255 for different number of samples
xi=10*cos(2*pi*ti*1e3)+6*cos(2*pi*ti*6e3);
k=4;
Fs=k*Fi;
t=[0:1/Fs:255/Fi+(k-1)/Fs];
L=length(t);
for i=1:1:length(ti)
    x(k*(i-1)+1)=xi(i);
    for j=1:1:k-1
        x(k*(i-1)+1+j)=0;
    end
end
Y=fft(x);
[B,A]=butter(5,1/k);
x=filter(B,A,x);
subplot(2,1,1);
plot(t,x);
title('Time Domain Analysis');
xlabel('time(secs)');
ylabel('Amplitute');

X=fft(x);
% F=[0:Fs/L:Fs/2-Fs/L,-Fs/2:Fs/L:-Fs/L];
%plot(F,2*abs(X)/L);
F=-Fs/2:Fs/L:Fs/2-Fs/L;
X=fftshift(X);    
subplot(2,2,3);
plot(F,2*abs(Y)/L);
title('Magnitute Spectrum (Before LPF)');
xlabel('Frequency(Hz)');
ylabel('Magnitute');

subplot(2,2,4);
plot(F,2*abs(X)/L);
title('Magnitute Spectrum (After LPF)');
xlabel('Frequency(Hz)');
ylabel('Magnitute');