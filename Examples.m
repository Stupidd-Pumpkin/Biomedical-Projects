% x=[0:0.1:1];
% y=x.*exp(-2.*x);
% for i=1:11
%     if y(i)==max(y)
%         result=(i-1)/10;
%     end
% end

% q2function [oddSum,funVal] =a1p2(x);
% oddSum=0;
% funVal=0;
% for i=1:2:numel(x)
%     oddSum=oddSum+x(i);
%     funVal=funVal+x(i);
% end
% for i=2:2:numel(x)
%     funVal=funVal-x(i);
% end

% % CurrentA=12;
% % CurrentB=15;
% % A=[CurrentA];
% % B=[CurrentB];
% % i=1;
% % while(CurrentA<CurrentB)
% %     i=i+1;
% %     CurrentA=CurrentA*(106/100);
% %     CurrentB=CurrentB*(102/100);
% %     A=[A,CurrentA];
% %     B=[B,CurrentB];
% % end

% % % t=[0:0.04:4];
% % % time=[0:0.4:4];
% % % y=20.*t-9.8.*(t.^2)/2;
% % % location=20.*time-9.8.*(time.^2)/2;
% % % plot(t,y,'-b');
% % % hold on;
% % % plot(time,location,'bo');



% x=1;
% truevalue=1/(1+x^2);
% for i=1:10
%     h=10^(-i);
%     hall(i)=h;
%     aprvalue(i)=(4*atan(x+h)-atan(x+2*h)-3*atan(x))/(2*h);
%     err(i)=abs(truevalue-aprvalue(i));
%     if err(i)==min(err)
%         hmin=h;
%         errmin=err(i);
%     end
% end
% loglog(hall,err);
%     hmin
%     errmin
    
% function [x,y,n] = myCubeRoot(a,tol)
% trueval=a^(1/3);
% x=1/2;
% n=0;
% y=x;
% tolerence=abs(trueval-x);
% while (tolerence>tol)
%     x=0.5*(x+(a/(x^2)));
%     tolerence=abs(trueval-x);
%     n=n+1;
%     y=[y,x];
% end

% % function cosval=a2p3(x,n)
% % vec=[1,(x.^(1:n))./cumprod(1:n),0,0];
% % i=1;
% % cosval=1;
% % while ((4*i)-2)<=n
% %     cosval=cosval-vec((4*i)-1)+vec((4*i)+1);
% %     i=i+1;
% % end

% % % disp('F off');



% a=0;
% b=2;
% n=20000000;
% x=[a:((b-a)/n):b];
% simpsons=((b-a)/(3*n))*(a3p1f(a)+a3p1f(b));
% for i=1:n
%     if mod(i,2)==1
%     simpsons=simpsons+((b-a)/(3*n))*(4*(a3p1f(x(i))));
%     else
%         simpsons=simpsons+((b-a)/(3*n))*(2*(a3p1f(x(i))));
%     end
% end

% function y=f(x)
% y=exp(x.^2);
% end
% 
% function y=a3p2(x)
% y=5*x.*log(x.^2);
% %%call using quad(@a3p2,0,2)
%
% function y=a3p3f(x)
% y=0.2*x.*log(x.^2);


% % a=2;
% % b=3;
% % h=0.0002;
% % x=[a:h:b];
% % trapezoid=(h/2)*[-a3p3f(a)+2*sum(a3p3f(x))-a3p3f(b)];




% % a=[0];
% % u=10;
% % t=0;
% % i=1;
% % g=9.8;
% % y=0;a
% % time=0;
% % location=0;
% % while(y>=0)
% %     a(i)=y;
% %     disp(['At t=',num2str(t),'      y=',num2str(y)]);
% %     t=t+0.1;
% %     y=u*t-g*t*t/2;
% %     i=i+1;
% %     time=[time;t];
% %     location=[location;y];
% % end
% % a
% % plot(time,location,'-.bx');
% % grid on;
    