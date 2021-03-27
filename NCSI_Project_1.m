clear all; close all;
fiberType = 3; % spontaneous rate of the fiber; "1" = Low; "2" = Medium; "3" = High
%F0=CF;
T  = 200e-3;  % stimulus duration in seconds
rt = 10e-3;   % rise/fall time in seconds
stimdb = 80; % stimulus intensity in dB SPL
nrep = 50;               % number of stimulus repetitions (e.g., 50)
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
t = 0:1/Fs:T-1/Fs; % time vector

%%
%%Part 1
CF = 500; % CF in Hz;
R_I1(1)=0;
for i=0:9
    clear R
    R(1)=0;
    stimdb = 10*(i-1);
    for j=0:1:72
        F0 = 62.5*2.^(j/8);     % stimulus frequency in Hz
        pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
        [vihc,synout,Psth,psthtime,timeout] = Model(CF,fiberType,T,rt,nrep,pin,t,Fs);
        R(j+1)=mean(Psth(1:floor(end/2)))-mean(Psth(ceil(end/2):end));
        if F0==500
            R_I1(i+1)=R(j+1);
        end
    end
    plot(62.5*2.^(0:1/8:9),R)
    hold on
    %     R(i+1)=numel(findpeaks(Psth(1:floor(end/2))))- numel(findpeaks(Psth(ceil(end/2):end)));
end
title('Frequency response of ANF with BF=500Hz');
xlabel('Frequency, Hz');ylabel('Average rate');
legend('-10dB', '0dB','10dB','20dB','30dB','40dB','50dB','60dB','70dB','80dB');
hold off

figure()
CF  = 4000; % CF in Hz;
R_I2(1)=0;
for i=0:9
    clear R
    R(1)=0;
    stimdb = 10*(i-1);
    for j=0:1:72
        F0 = 62.5*2.^(j/8);     % stimulus frequency in Hz
        t = 0:1/Fs:T-1/Fs; % time vector
        pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
        [vihc,synout,Psth,psthtime,timeout] = Model(CF,fiberType,T,rt,nrep,pin,t,Fs);
        R(j+1)=mean(Psth(1:floor(end/2)))-mean(Psth(ceil(end/2):end));
        if F0==4000
            R_I2(i+1)=R(j+1);
        end
    end
    plot(62.5*2.^(0:1/8:9),R)
    hold on
    %     R(i+1)=numel(findpeaks(Psth(1:floor(end/2))))- numel(findpeaks(Psth(ceil(end/2):end)));
end
title('Frequency response of ANF with BF=4000Hz');
xlabel('Frequency, Hz');ylabel('Average rate')
legend('-10dB', '0dB','10dB','20dB','30dB','40dB','50dB','60dB','70dB','80dB');
hold off

figure()
plot(-10:10:80,R_I1)
hold on
plot(-10:10:80,R_I2)
title('Rate vs Intensity');
xlabel('Intensity (dB SPL)');ylabel('Average rate at Best frequency')
legend('BF=500', 'BF=4000');

%%
%%part 2
clear pin
clear R
[fivewo,Fs]=audioread('fivewo.wav');
fivewo=resample(fivewo,1024,1000);
fivewo=fivewo';
pin=fivewo(105001:115500);
Fs=Fs*1024/1000;
% sound(pin,Fs);
CF=500;
T=length(pin)/Fs;
rt = T/10;
t = 0:1/Fs:T-1/Fs;
R(1)=0;
for i=1:11
    p=10^((i-13)*0.5)*2*pin/rms(pin);
    stimdb=20*log10(rms(p)*(10^6)/20);
    [vihc,synout,Psth,psthtime,timeout] = Model(CF,fiberType,T,rt,nrep,p,t,Fs);
    R(i)=mean(Psth(1:floor(end/2)))-mean(Psth(ceil(end/2):end));
end

figure()
plot(-20:10:80, R)
title('Rate Intensity function for the vowel "ah"');
xlabel('Intensity (dB SPL)');ylabel('Average rate response')
% legend('');


clear R
pin=fivewo(1:157500);
T=length(pin)/Fs;
rt = T/10;
t = 0:1/Fs:T-1/Fs;
R(1,1,1)=0;
Ravg(1,1)=0;
i=0;
for q=[6,7,9] %[30,40,60] dB spl
    p=10^((q-13)*0.5)*2*pin/rms(pin);
    stimdb=20*log10(rms(p)*(10^6)/20);
    i=i+1;
    for j=1:55
        CF = 62.5*2.^(j/8+0.25);
        [vihc,synout,Psth,psthtime,timeout] = Model(CF,fiberType,T,rt,nrep,p,t,Fs);
        for k=1:floor(length(Psth)/2)
            R(i,j,k)=Psth(k);
        end
        Ravg(i,j)=mean(Psth(1:floor(end/2)))-mean(Psth(ceil(end/2):end));
    end
    
end

figure()
hold on
plot(62.5*2.^(3/8:1/8:57/8),Ravg(1,:));
plot(62.5*2.^(3/8:1/8:57/8),Ravg(2,:));
plot(62.5*2.^(3/8:1/8:57/8),Ravg(3,:));
% title('Histogram');
xlabel('Frequency (Hz)');ylabel('Average rate')
legend('Intensity=30dB','Intensity=40dB','Intensity=60dB');
title('Rate response to the speech');
hold off


figure()
spectrogram(pin,hanning(floor(25.6*Fs/1000)),floor(12.8*Fs/1000),1024,Fs,'yaxis')
title('Spectrogram of the input signal');

for q2=1:3
    for q=2:7
        w=2*2^q;
        clear Rwavg
        Rwavg(1,1)=0;
        for k=1:55
            j=0;
            for i=0:w/2:length(R(1,1,:))-w
                j=j+1;
                Rwavg(k,j)=mean(R(q2,k,i+1:i+w));
            end
        end
        figure()
        imagesc([0 T],[1 55],Rwavg);
        ticks=62.5*2.^((1:5:55)/8+0.25);
        set(gca, 'YDir', 'normal','YTick',1:5:55,'YTickLabel',ticks);
        xlabel('Time (secs)');ylabel('Frequency (Hz) -log scale')
        title(['Rate response vs Best frequencies with a window of ' num2str(w/2) ' ms'])
        colorbar
    end
    
end
%%
%%part 3

% mxpts = length(t);
% irpts = round(rt*Fs);
% pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
% pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;


for q2=[6,7,9]
    clear PL_F
    clear PL_V
    PL_F(1,1)=0;
    PL_V(1,1)=0;
    for k=0:8
        CF=125*2^(k/2);
        p=10^((q2-13)*0.5)*2*pin/rms(pin);
        vihc = catmodel_IHC(p,CF,nrep,1/Fs,T,1,1);
        [synout,Psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,0);
        psthtime = (1:length(Psth))*1/Fs;
        w=100*12.8;
        j=0;
        for i=0:w/2:length(Psth)-w
            j=j+1;
            PSTH=fft(Psth(i+1:i+w));
            L=length(PSTH);
            PSTH = abs(PSTH(1:L/2+1)/L);
            PSTH(2:end-1) = 2*PSTH(2:end-1);
            f = Fs*(0:L/2)/L;
            [PL_V(k+1,j),q]=max(PSTH(2:end));
            PL_F(k+1,j)=f(q(1)+1);
        end
    end
    figure()
    spectrogram(pin,hanning(floor(25.6*Fs/1000)),floor(12.8*Fs/1000),1024,Fs,'yaxis')
    colorbar
    hold on
    for i=1:9
        plot(length(t)/(2*Fs*length(PL_F(1,:))):length(t)/(Fs*length(PL_F(1,:))):length(t)/Fs,PL_F(i,:)/1000,'*')
    end
    hold off
    legend('Spectrogram','BF=125 Hz','BF=176.8 Hz','BF=250 Hz','BF=353.6 Hz','BF=500 Hz','BF=707.1 Hz','BF=1000 Hz','BF=1414.2 Hz','BF=2000 Hz')
    
    figure()
    hold on
    for i=1:9
        plot(PL_V(i,:))
    end
    hold off
end
