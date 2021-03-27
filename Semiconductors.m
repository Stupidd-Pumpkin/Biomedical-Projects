%%Roll_NO :15EC10054
%Defining Constants
h = 6.625e-34;
m = 9.11e-31;
q = 1.602e-19;
a = 14e-9;
b1 = 15e-9;
b2 = 7.5e-9;
d = 0.021e-9;
n=2;
C=(h^2)/(8*pi*pi*m);
b=b1;
%b=b2 for the second case
N=floor(((n+1)*a+n*b)/d)+1;
prin_diag=2*C/(d^2);
side_diag=-1*C/(d^2);
x=0;
U0=1e-20;
for i=1:N-2
    x=[x,x(end)+d];
    for j=0:n
        if (x(end) > (j*a+(j-1)*b) && x(end) < (j*(a+b)))
            potential = U0;
        else
            potential=0;
        end
    end
    prin_diag=[prin_diag,2*C/(d^2)+potential];
    side_diag=[side_diag,-C/(d^2)];
end
 x=[x,x(end)+d];
 prin_diag=[prin_diag,2*C/(d^2)];
H=diag(prin_diag)+diag(side_diag,-1)+diag(side_diag,1);
[V,E]=eig(H);
%Plotting The Wave Function with respect to distance
for i=1:4
plot(x,V(1:N,i));
hold on;
end
hold off
grid on



% %%Roll_NO :15EC10054
% %Defining Constants
% h = 6.625e-34;
% m = 9.11e-31;
% q = 1.602e-19;
% a = 14e-9;
% b1 = 15e-9;
% b2 = 7.5e-9;
% d = 0.021e-9;
% n=1;
% C=(h^2)/(8*pi*pi*m);
% b=b1;
% %b=b2 for the second case
% N=floor(((n+1)*a+n*b)/d)+1;
% prin_diag=2*C/(d^2);
% side_diag=-1*C/(d^2);
% x=0;
% U0=1e-20;
% for i=1:N-2
%     x=[x,x(end)+d];
%     for j=1:n
%         if (x(end) > (j*a+(j-1)*b) && x(end) < (j*(a+b)))
%             potential = U0;
%         else
%             potential=0;
%         end
%     end
%     prin_diag=[prin_diag,2*C/(d^2)+potential];
%     side_diag=[side_diag,-C/(d^2)];
% end
%  x=[x,x(end)+d];
%  prin_diag=[prin_diag,2*C/(d^2)];
% H=diag(prin_diag)+diag(side_diag,-1)+diag(side_diag,1);
% [V,E]=eig(H);
% %Plotting The Wave Function with respect to distance
% for i=1:100
% plot(x,V(1:N,i));
% hold on;
% end
% hold off
% grid on
% 



% % h = 6.625e-34;
% % m = 9.11e-31;
% % q = 1.602e-19;
% % L=1e-18;
% % k=h^2/(8*(pi^2)*m);
% % N=9;
% % x=[-0.4:0.1:0.4];
% % d=L/(2*N);
% % prin_diag=2*k/(d^2);
% % side_diag=[-k/(d^2)+0*[1:8]];
% % prin_diag=[2*k/(d^2)+0.1*[4:-1:1],prin_diag,2*k/(d^2)+0.1*[1:1:4]];
% % H=diag(prin_diag)+diag(side_diag,-1)+diag(side_diag,1);
% % [V,E]=eig(H);
% % plot(x,V(1:N,1));

