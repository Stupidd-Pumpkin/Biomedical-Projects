function [T]= Bin2Ther(B,nob)
T=zeros(2^nob,1);
T(1:B)=1;