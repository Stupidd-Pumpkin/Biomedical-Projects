clear all;
Fs=8e3;
in=input('Give the input inside ' ' (as string):');
Ts=1000/Fs;
ts=(1/Fs:1/Fs:Ts);
t=0;
x=0;
for i=1:length(in)
    t=[t,ts+max(t)];
    if in(i)=='1'
        x=[x,cos(2*pi*697*ts)+cos(2*pi*1209*ts)];
    elseif in(i)=='2'
        x=[x,cos(2*pi*697*ts)+cos(2*pi*1336*ts)];
    elseif in(i)=='3'
        x=[x,cos(2*pi*697*ts)+cos(2*pi*1477*ts)];
    elseif in(i)=='A'
        x=[x,cos(2*pi*697*ts)+cos(2*pi*1633*ts)];
    elseif in(i)=='4'
        x=[x,cos(2*pi*770*ts)+cos(2*pi*1209*ts)];
    elseif in(i)=='5'
        x=[x,cos(2*pi*770*ts)+cos(2*pi*1336*ts)];
    elseif in(i)=='6'
        x=[x,cos(2*pi*770*ts)+cos(2*pi*1477*ts)];
    elseif in(i)=='B'
        x=[x,cos(2*pi*770*ts)+cos(2*pi*1633*ts)];
    elseif in(i)=='7'
        x=[x,cos(2*pi*852*ts)+cos(2*pi*1209*ts)];
    elseif in(i)=='8'
        x=[x,cos(2*pi*852*ts)+cos(2*pi*1336*ts)];
    elseif in(i)=='9'
        x=[x,cos(2*pi*852*ts)+cos(2*pi*1477*ts)];
    elseif in(i)=='C'
        x=[x,cos(2*pi*852*ts)+cos(2*pi*1633*ts)];
    elseif in(i)=='.'
        x=[x,cos(2*pi*941*ts)+cos(2*pi*1209*ts)];
    elseif in(i)=='0'
        x=[x,cos(2*pi*941*ts)+cos(2*pi*1336*ts)];
    elseif in(i)=='#'
        x=[x,cos(2*pi*941*ts)+cos(2*pi*1477*ts)];
    elseif in(i)=='D'
        x=[x,cos(2*pi*941*ts)+cos(2*pi*1633*ts)];
    end
end


x(1)=[];
t(1)=[];

L=length(t);
Ls=length(ts);
X=fftshift(fft(x));
F=-Fs/2:Fs/L:Fs/2-Fs/L;
f=-Fs/2:Fs/Ls:Fs/2-Fs/Ls;

x=awgn(x,0);
X=fftshift(fft(x));
F=-Fs/2:Fs/L:Fs/2-Fs/L;
f=-Fs/2:Fs/Ls:Fs/2-Fs/Ls;

xlswrite('DTMF.xlsx',x','x');
[q, qq, Raw]=xlsread('DTMF.xlsx','x');
[Raw{:, :}]=deal(NaN);
xlswrite('DTMF.xlsx', Raw,'x');
xlswrite('DTMF.xlsx', Raw,'t');

xlswrite('DTMF.xlsx',x','x');
xlswrite('DTMF.xlsx',t','t');


figure(1)
subplot(2,1,1);
plot(t,x);
title('Input Signal');
xlabel('time(secs)');
ylabel('Amplitute');
subplot(2,1,2);
plot(F,2*abs(X)/L);
title('Magnitute Spectrum');
xlabel('Frequency(Hz)');
ylabel('Magnitute');