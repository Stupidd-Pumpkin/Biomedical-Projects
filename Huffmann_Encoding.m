clear all;
x = rgb2gray(imread('003.png'));
s=size(x);
p(1:1:256)=0;

for i=1:1:s(1)
    for j=1:1:s(2)
        p(x(i,j)+1)=p(x(i,j)+1)+1;
        end
end
y(:,1)=1:1:256;
y(:,2)=p/(s(1)*s(2));
y(:,3)=0;y(:,4)=0;y(:,5)=0;

[~,index]=sort(y(:,2));
z=y(index,:);
%1 has 8 bit value, 2 has probabilities, 3 has left leaf, 4 has right leaf, 5 has bit assigned
for i=1:2:255*2-1
    node(2)=z(i,2)+z(i+1,2);
    node(3)=i;
    node(4)=i+1;
    node(1)=0;
    node(5)=0;
    for j=length(z):-1:1
        if node(2)>=z(j,2)
            z=[z(1:j,:);node(1,:);z(j+1:length(z),:)];
            break;
        end
    end
end

z=PostOrder(z,511,1);
j=1;
for i=1:1:511
    if z(i,3)==0
        f(j,1)=z(i,1);
        f(j,2)=z(i,2);
        f(j,3)=z(i,5)-128;
        j=j+1;
    end
end

q=sum(f(:,2).*log(f(:,3)/log(2)));
r=sum(f(:,2).*log(f(:,1)))/log(2);
% function matrix= PostOrder(matrix,index,x)
% if matrix(index,3)
%     matrix=PostOrder(matrix,matrix(index,3),2*x);
%     matrix=PostOrder(matrix,matrix(index,4),2*x+1);
% end
% matrix(index,5)=x;
% end

