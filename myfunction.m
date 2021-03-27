function [y] =Polynimial(n,x);
c=1./[1:n];
y=[1,c.*x.^c];
end
