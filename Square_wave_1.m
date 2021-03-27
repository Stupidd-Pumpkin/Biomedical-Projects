%Square wave
Fs=20; %KHz
t=0:1/Fs:255/Fs; %ms
L=length(t);
i=1;
x=0;
while(i<=length(t))
    if mod(ceil(t(i)/0.5),2)==1
        x(i)=1;
    else
        x(i)=-1;
    end
    i=i+1;
end
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