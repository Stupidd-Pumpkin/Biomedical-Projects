clear all;
image =(rgb2gray(imread('Team_7.png')));
x=double(image);
s=size(x);
B=2;

a(1:1:2^B)=2^(8-B):2^(8-B):2^8;
b(1)=0;
b(2^B+1)=255;
for i=1:1:length(a)-1
    b(i+1)=(a(i)+a(i+1))/2;
end

for k=1:1:10
    for r=1:1:2^B
        sum_a=0;
        count=0;
        for i=1:1:s(1)
            for j=1:1:s(2)
                if x(i,j)>=b(r) && x(i,j)<=b(r+1)
                    sum_a=sum_a+x(i,j);
                    count=count+1;
                end
            end
        end
        a(r)=floor(sum_a/count);
        if count==0
            a(r)=a(r-1);
        end
    end
    b(1)=0;
    for i=1:1:length(a)-1
        b(i+1)=(a(i)+a(i+1))/2;
    end
    b(2^B+1)=255;
end

a=floor(a);
b=floor(b);
MSE=0;

for i=1:1:s(1)
    error=0;
    for j=1:1:s(2)
        if x(i,j)+error<=b(1)
            y(i,j)=a(1);
            MSE=MSE+(y(i,j)-x(i,j))^2;
        end
        if x(i,j)+error>=b(2^B+1)
            y(i,j)=a(2^B);
            MSE=MSE+(y(i,j)-x(i,j))^2;
        end
        for r=1:1:2^B
            if x(i,j)+error>=b(r) && x(i,j)+error<b(r+1)
                y(i,j)=a(r);
                MSE=MSE+(y(i,j)-x(i,j))^2;
            end
        end
        error=x(i,j)+error-y(i,j);
        if abs(error)>=255
            error=0;
        end
    end
end
y=uint8(y);
imwrite(y,'DFTBA_Super_Cool.jpg');