% amp=[];
% phase=[0:pi/50:10*pi];
% x=0;
% for x=0:pi/50:10*pi
%     Y=0;
%     for i=1:100
%         y(i)=sin((2*i-1)*x)/(2*i-1);
%         Y=Y+y(i);
%     end
%     amp=[amp,Y];
% end
% plot(phase,amp);

amp=[];
time=[0:100];
t=0;
for t=0:100
    Y=0;
    for i=1:1000
        y(i)=sin((2*i-1)*(t+10))/(2*i-1);
        Y=Y+y(i);
    end
    amp=[amp,Y];
end
plot(time,amp);