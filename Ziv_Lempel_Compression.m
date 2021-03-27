clear all;
text=load('experimental.txt');


B=8;
p=2^B-1;
%%%%%%%%%%%%%%%
%Forming Dictionary
dict(1:p,1:3)=0;
dict(1:p,1)=1:p;
k=0;
for i=1:p
    k=k+1;
    dict(i,2)=2+data(k);
    dict(i,3)=data(k);
    for j=1:i-1
        if dict(i,2)==dict(j,2)
            k=k+1;
            dict(i,2)=dict(i,2)*2+data(k);
            dict(i,3)=dict(j,1)*2+data(k);
        end
    end
end


%%%%%%%%%%%%%%%
%Encoding
out=dict(:,3);
i=p;
while(k<length(data))
    i=i+1;
    k=k+1;
    dat=2+data(k);
    for j=1:p
        if dat==dict(j,2)
            if k==length(data)
                break
            end
            k=k+1;
            dat=dat*2+data(k);
            flag=j;
        end
    end
    if k==length(data)
        break
    end
    out(i)=dict(flag,1)*2+data(k);
end

out=dec2bin(out);
x='';
for i=1:length(out)
    x=strcat(x,out(i,:));
    
end
y=double(x)-48;


%  file=fopen('Experimental_2.txt','w');
%  for i=1:length(y)
%      fprintf(file, '%f ',y(i));
%  end
%  fclose(file);

% for i=1:p
%     dict(i,2)=dict(i,2)-2^(floor(log(dict(i,2))/log(2)));
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Receiving End
B=8;
len=length(out);
dict_2(1:len,1)=1:len;
dict_2(1:len,3)=1;
for i=1:len
    q=bin2dec(string(y((B+1)*(i-1)+1:(B+1)*(i-1)+B+1)));
    for j=1:B+1
        dict_2(i,3)=dict_2(i,3)*2+q(j);
    end
end
dict_2(:,3)=dict_2(:,3)-512;

for i=1:len
    if floor(dict_2(i,3)/2)==0
        dict_2(i,2)=2+mod(dict_2(i,3),2);
    else
        dict_2(i,2)=dict_2(floor(dict_2(i,3)/2),2)*2+mod(dict_2(i,3),2);
    end
end

for i=1:length(out)
    if floor(dict_2(i,3)/2)==0
        dict_2(i,3)=2+mod(dict_2(i,3),2);
    else
        dict_2(i,2)=dict_2(floor(dict_2(i,3)/2),2)*2+mod(dict_2(i,3),2);
    end
end

DFTBA=de2bi(dict_2(:,2));
s=size(DFTBA);
l=0;
for i=1:s(1)
    j=s(2);
    while(~DFTBA(i,j))
        j=j-1;
    end
    for k=j-1:-1:1
        l=l+1;
        data(2,l)=DFTBA(i,k);
    end
end

error=0;
for i=1:length(data)
    if data(1,i)~=data(2,i)
        error=error+1;
    end
end

eff=(length(data)-length(y))/length(data);