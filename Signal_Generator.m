clear all

F=0.3e6; %frequency of the signal
Fs=F*100; %sampling frequency
N=5; %number of cycles
A=0.5; %amplitude
phi=-90; %phase in degrees
DC=0; %DC offset

t = 0:1/Fs:N/F-1/Fs;
L=length(t);
f=-Fs/2:Fs/L:Fs/2-Fs/L;

%sinewave
x=A*sin(2*pi*F*t+phi*pi/180)+DC;

%square wave
d= 50; %duty cycle
x=A*square(2*pi*F*t+phi*pi/180,d)+DC;

% % %Hamming Wave
% h=sin(2*pi*F/2*t+phi*pi/180);
% w=0.54-0.46*cos(2*pi*F*t/N);
% x=A*h.*w+DC;
% 
% [B,A]=butter(6,2*F/Fs);
% y=filter(B,A,x);

X=fftshift(fft(x));
figure(1)
% subplot(2,1,1)
plot(t,x,'LineWidth',5)
ylim([-0.7 0.7])
% subplot(2,1,2)
% plot(f,abs(X)*2/L)
