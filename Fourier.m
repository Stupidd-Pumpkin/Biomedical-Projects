x=-6*pi:pi/100:6*pi;
y=0;
for i=1:length(x)
    y(i)=0;
    for n=0.01:0.01:100
        y(i)=y(i)+exp(j*n*x(i))/(j*n);
    end
end
% plot(x,y,'k');
(max(y))
z=0;
for i=1:length(x)
    z(i)=0;
    for n=0.01:0.01:100
        z(i)=z(i)-exp(j*n*(x(i)-5))/(j*n);
    end
end
% hold on;
plot(x,y+z,'r');
