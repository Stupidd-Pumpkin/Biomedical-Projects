%Sine wave
f=1.5
Fs=50;
t=0:1/Fs:1000/Fs;
L=length(t);
i=1;
while(i<=L)
    x(i)=sin(2*pi*f*t(i));
    i=i+1;
end
plot(t,x);
X=fft(x);
F=[0:Fs/L:Fs/2-Fs/L,-Fs/2:Fs/L:0];
plot(F,2*abs(X)/L);