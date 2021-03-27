%Audio_Producer
Fs=12*10^3;
t=[0:1/Fs:10000*1/Fs]; %63,127,255 for different number of samples
L=length(t);
x=10*cos(2*pi*t*2*10^3);

y=audioplayer(x,Fs);
play(y);
plot(t,x);