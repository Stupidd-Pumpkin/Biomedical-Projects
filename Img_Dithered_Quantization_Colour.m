clear all;
img=imread('Team_7.PNG');
s=size(img);
for k=1:1:s(3)
    for i=1:1:s(1)
        error=0;
        for j=1:1:s(2)
            if img(i,j,k)+error>127
                y(i,j,k)=255;
            else
           y(i,j,k)=0;
       end
       error=(img(i,j,k)+error)-y(i,j,k);
    end
end
end
imwrite(y,'DFTBA_Dithered_Quantization_Colour.jpg');

% x=img(1:s(1),1:s(2),1);
% img(:,:,3)=img(:,:,1);
% img(:,:,1)=x;
