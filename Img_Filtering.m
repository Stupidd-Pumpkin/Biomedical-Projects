clear all;
x = rgb2gray(imread('Apple.JPG'));
X=fftshift(fft2(double(x)));
X2=mat2gray(log(abs(X)+1));
figure(1);
imshow(X2);
s=size(X);
Fc=min(s(1),s(2))/16;  %the constant can be changed to shift the cut off frequency of the filter used
Y=X;
for i=1:s(1)
    for j=1:s(2)
        if (s(1)/2-i)^2+(s(2)/2-j)^2<=Fc^2  
            % <= behaves like HPF, and >= behaves like LPF
            Y(i,j)=0;
        end
    end
end
figure(2)
% imshow(Y);
y=ifft2(fftshift(Y));
% y=uint8(abs(y));
imshow(y);

kernel=[-1,-1,-1;-1,8,-1;,-1,-1,-1];
Z=uint8(conv2(x,kernel,'same'));
figure(3);
imshow(Z);
