%Harmonics and MSE 
clear all;
Fs=1000; %KHz
t=0:1/Fs:10; %ms
L=length(t);
i=1;
while(i<=length(t))
    if mod(ceil(t(i) /0.5),2)==1
        x(i)=1;
    else
        x(i)=-1;
    end
    i=i+1;
end
Fsi=10;nhhu
k=floor(Fs/Fsi);
t1=[0:1/Fsi:floor(length(t)/k)/Fsi];
for i=1:1:length(t1)
    x1(i)=x(k*i-(k-1));
end

for i=1:1:length(t1)
    x2(k*(i-1)+1)=x1(i);
end

[B,A]=butter(5,1/k);
x2=k*filter(B,A,x2);
X=fft(x2);
F=-Fs/2:Fs/L:Fs/2-Fs/L;
Y=fftshift(X);    
y=ifft(X);
z=(x-x2).^2/L;
s=sum(z);


figure(1);
subplot(2,1,1);
plot(t,x);
title('Original Signal');
xlabel('time(msecs)');
ylabel('Amplitute');
subplot(2,1,2);
plot(t,x2);
title('Signal after Sampling');
xlabel('time(msecs)');
ylabel('Amplitude');

figure(2);
subplot(2,1,1);
plot(F,2*abs(Y)/L);
title('Magnitute Spectrum');
xlabel('Frequency(KHz)');
ylabel('Magnitute');
subplot(2,1,2);
plot(t,z);
title('Error');
xlabel('time(msecs)');
ylabel('Amplitude');


