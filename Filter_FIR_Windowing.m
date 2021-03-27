clear all;
wc1=0.05*pi;
wc2=0*pi;
N=8;
k=(N-1)/2;
n=0:N-1;
hd1(1:N)=0;
hd2(1:N)=0;
Fi=-pi:2*pi/N:pi-2*pi/N;

for i=0:N-1
    if i==k
        hd1(i+1)=wc1/pi;
    else
        hd1(i+1)=sin(wc1*(i-k))/(pi*(i-k));
    end
end

for i=0:N-1
    if i==k
        hd2(i+1)=wc2/pi;
    else
        hd2(i+1)=sin(wc2*(i-k))/(pi*(i-k));
    end
end


% % [Hd,F]=freqz(hd,1,Fi);
% Hd=fftshift(fft(hd1));
% figure(1);
% subplot(2,1,1);
% plot(n,hd1);
% title('hd');
% subplot(2,1,2);
% plot(Fi,abs(Hd));

%%Rectangular window
w(n+1)=1;

%Triangular (Bartlett) Window
w(n+1)=(n/k).*(n<=k)+(2-n/k).*(n>k);

%Hanning Window
w(n+1)=0.5-0.5*cos(pi*n/k);

%Hamming Window
w(n+1)=0.54-0.46*cos(pi*n/k);
% 
% %blackman Window
% w(n+1)=0.42-0.5*cos(pi*n/k)+0.08*cos(2*pi*n/k);
% 
% %Kaiser Window
% beta=1;
% w(n+1)=kaiser(N,beta);

figure(2);
title('Window function');
freqz(w,1,Fi);
% subplot(2,1,1);
plot(n,w/max(w),'LineWidth',4);
hold on
plot(n,ones(8,1),'LineWidth',4);
plot(n,1/10*[4,6,8,10,10,8,6,4],'LineWidth',4);
ylim([0 1.2])
hold off
% subplot(2,1,2);
% W=fftshift(fft(w));
% plot(Fi,mag2db(abs(W)));
% 
% h=(hd1-hd2).*w;
% figure(3)
% freqz(h);
% 
% figure(4)
% subplot(2,1,1);
% plot(n,h);
% subplot(2,1,2);
% H=fftshift(fft(h));
% plot(Fi,abs(H));

% Fs=100e3;
% t=(0:1/Fs:2048/Fs); %63,127,255 for different number of samples
% L=length(t);
% F=-Fs/2:Fs/L:Fs/2-Fs/L;
% F=2*pi*F/Fs;
% x=10*cos(2*pi*t*2e3)+6*cos(2*pi*t*3e3)+8*cos(2*pi*t*4e3);
% X=freqz(x,1,F);
% X=fftshift(fft(x));
% 
% figure(5);
% subplot(2,1,1);
% plot(t,x);
% title('Input Signal');
% subplot(2,1,2);
% plot(F,mag2db(abs(X)));
% 
% noise=10*rand([1,L])-5;
% NOISE=freqz(noise,1,F);
% 
% figure(6);
% subplot(2,1,1);
% plot(t,x+noise);
% title('Noisy Input Signal');
% subplot(2,1,2);
% plot(F,mag2db(abs(X+NOISE)));
% 
% y=filtfilt(h,1,x);
% figure(7)
% freqz(y,1,F);
% 
% figure(8);
% subplot(2,1,1);
% plot(t,y);
% title('Output Signal');
% subplot(2,1,2);
% plot(F,(abs(fftshift(fft(y)))));
