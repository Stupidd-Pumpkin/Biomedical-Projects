clear all;
x = rgb2gray(imread('Team_7.png'));
s=size(x);
bits=4;
p=2^bits;
y=0;
for i=1:1:s(1)
    for j=1:1:s(2)
        for k=1:1:p
            if ((x(i,j)>=256*(k-1)/p)&&(x(i,j)<256*(k)/p))
                y(i,j)=255*(k-1)/p;
            end
       end
    end
end
y=uint8(y);
imwrite(y,'DFTBA3.jpg')
info=imfinfo('DFTBA3.jpg');
imshow(y);
xd = whos('x');
yd = whos('y');
fprintf('grey: %f\n b/w: %f\n', xd.bytes, yd.bytes);
