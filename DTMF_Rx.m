clear all;
Fs=8e3;
x=xlsread('DTMF.xlsx','x');
t=xlsread('DTMF.xlsx','t');

Ts=1000/Fs;
ts=(1/Fs:1/Fs:Ts);
L=length(t);
Ls=length(ts);

X=fftshift(fft(x));
F=-Fs/2:Fs/L:Fs/2-Fs/L;
f=-Fs/2:Fs/Ls:Fs/2-Fs/Ls;

[N,Fo,Ao,W]=firpmord([1,1.1],[1,0],[0.1,0.1],8);
h=firpm(N,Fo,Ao,W);
y1=filtfilt(h,1,x);
Y1=fftshift(fft(y1));

[N,Fo,Ao,W]=firpmord([1,1.1],[0,1],[0.1,0.1],8);
h=firpm(N,Fo,Ao,W);
y2=filtfilt(h,1,x);
Y2=fftshift(fft(y2));

%%%%%%%%%%%%%%%
[N,Fo,Ao,W]=firpmord([0.65,0.68,0.72,0.75],[0,1,0],[0.1,0.1,0.1],8);
h=firpm(N,Fo,Ao,W);
y1_1=filtfilt(h,1,y1);
Y1_1=fftshift(fft(y1_1));

[N,Fo,Ao,W]=firpmord([0.72,0.75,0.79,0.82],[0,1,0],[0.1,0.1,0.1],8);
h=firpm(N,Fo,Ao,W);
y1_2=filtfilt(h,1,y1);
Y1_2=fftshift(fft(y1_2));

[N,Fo,Ao,W]=firpmord([0.8,0.83,0.87,0.9],[0,1,0],[0.1,0.1,0.1],8);
h=firpm(N,Fo,Ao,W);
y1_3=filtfilt(h,1,y1);
Y1_3=fftshift(fft(y1_3));

[N,Fo,Ao,W]=firpmord([0.89,0.92,0.96,0.99],[0,1,0],[0.1,0.1,0.1],8);
h=firpm(N,Fo,Ao,W);
y1_4=filtfilt(h,1,y1);
Y1_4=fftshift(fft(y1_4));

%%%%%%%%%%%%%%%%
[N,Fo,Ao,W]=firpmord([1.12,1.17,1.23,1.28],[0,1,0],[0.2,0.2,0.2],8);
h=firpm(N,Fo,Ao,W);
y2_1=filtfilt(h,1,y2);
Y2_1=fftshift(fft(y2_1));

[N,Fo,Ao,W]=firpmord([1.25,1.3,1.36,1.41],[0,1,0],[0.2,0.2,0.2],8);
h=firpm(N,Fo,Ao,W);
y2_2=filtfilt(h,1,y2);
Y2_2=fftshift(fft(y2_2));

[N,Fo,Ao,W]=firpmord([1.4,1.45,1.51,1.56],[0,1,0],[0.2,0.2,0.2],8);
h=firpm(N,Fo,Ao,W);
y2_3=filtfilt(h,1,y2);
Y2_3=fftshift(fft(y2_3));

[N,Fo,Ao,W]=firpmord([1.55,1.6,1.66,1.71],[0,1,0],[0.2,0.2,0.2],8);
h=firpm(N,Fo,Ao,W);
y2_4=filtfilt(h,1,y2);
Y2_4=fftshift(fft(y2_4));

%%%%%%%%%%%%%%%%


for i=1:L/Ls
    E1(i)=sum((y1_1((i-1)*Ls+1:Ls*i)).^2);
    E2(i)=sum((y1_2((i-1)*Ls+1:Ls*i)).^2);
    E3(i)=sum((y1_3((i-1)*Ls+1:Ls*i)).^2);
    E4(i)=sum((y1_4((i-1)*Ls+1:Ls*i)).^2);
    if E1(i)>E2(i) && E1(i)>E3(i) && E1(i)>E4(i)
        logic(i)=0;
    elseif E2(i)>E3(i) && E2(i)>E4(i)
        logic(i)=4;
    elseif E3(i)>E4(i)
        logic(i)=8;
    else
        logic(i)=12;
    end
end

for i=1:L/Ls
    E1(i)=sum((y2_1((i-1)*Ls+1:Ls*i)).^2);
    E2(i)=sum((y2_2((i-1)*Ls+1:Ls*i)).^2);
    E3(i)=sum((y2_3((i-1)*Ls+1:Ls*i)).^2);
    E4(i)=sum((y2_4((i-1)*Ls+1:Ls*i)).^2);
    if E1(i)>E2(i) && E1(i)>E3(i) && E1(i)>E4(i)
        logic(i)=logic(i)+0;
    elseif E2(i)>E3(i) && E2(i)>E4(i)
        logic(i)=logic(i)+1;
    elseif E3(i)>E4(i)
        logic(i)=logic(i)+2;
    else
        logic(i)=logic(i)+3;
    end
end
z='';
for i=1:L/Ls
    if logic(i)==0
        z=[z,'1'];
    elseif logic(i)==1
        z=[z,'2'];
    elseif logic(i)==2
        z=[z,'3'];
    elseif logic(i)==3
        z=[z,'A'];
    elseif logic(i)==4
        z=[z,'4'];
    elseif logic(i)==5
        z=[z,'5'];
    elseif logic(i)==6
        z=[z,'6'];
    elseif logic(i)==7
        z=[z,'B'];
    elseif logic(i)==8
        z=[z,'7'];
    elseif logic(i)==9
        z=[z,'8'];
    elseif logic(i)==10
        z=[z,'9'];
    elseif logic(i)==11
        z=[z,'C'];
    elseif logic(i)==12
        z=[z,'.'];
    elseif logic(i)==13
        z=[z,'0'];
    elseif logic(i)==14
        z=[z,'#'];
    elseif logic(i)==15
        z=[z,'D'];
    end
end

fprintf('The input is:\n');
disp(z);