clear all;
img=rgb2gray(imread('Team_7.PNG'));
s=size(img);
for i=1:1:s(1)
    error=0;
    for j=1:1:s(2)
       if img(i,j)+error>127
           y(i,j)=255;
       else
           y(i,j)=0;
       end
       error=(img(i,j)+error)-y(i,j);
    end
end
imwrite(y,'DFTBA_Dithered_Quantization.jpg');

% x=img(1:s(1),1:s(2),1);
% img(:,:,3)=img(:,:,1);
% img(:,:,1)=x;
