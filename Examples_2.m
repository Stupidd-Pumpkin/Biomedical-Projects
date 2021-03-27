f0=1e5;
w0=2*pi*f0;
f=[10:100:2e5];
w=2*pi*f;
V=50;
Rs=50;
Rl=1000;
L=((Rl-Rs)*Rs)^0.5/w0;
C=(Rl-Rs)^0.5/(Rl*Rs^0.5*w0);
Xl=w*L*1j;
Xc=-1j./(w*C);
Z=Rs+Xl+((Xc*Rl)./(Rl+Xc));
I=V./Z;
Io=I.*Xc./(Rl+Xc);
P=abs(Io.^2)*Rl;
plot(f,P,'-r');
Phalf=max(P)/2;
BW=0;
for i=1:2000
    if P(i) > Phalf
        BW=[BW,f(i)];
    end
end
bandwidth=BW(end)-BW(2);


% f0=1e5;
% w0=2*pi*f0;
% f=[10:100:2e5];
% w=2*pi*f;
% V=50;
% Rs=50;
% Rl=1000;
% C = 1/(w0*(sqrt(Rs*(Rl-Rs))));
% L = (Rl/w0)*(sqrt(Rs/(Rl-Rs)));
% Xl=w*L*1j;
% Xc=-1j./(w*C);
% Z=Rs+Xc+((Xl*Rl)./(Rl+Xl));
% I=V./Z;
% Io=I.*Xl./(Rl+Xl);
% P=abs(Io.^2)*Rl;
% plot(f,P,'-b');
% Phalf=max(P)/2;
% BW=0;
% for i=1:2000
%     if P(i) > Phalf
%         BW=[BW,f(i)];
%     end
% end
% bandwidth=BW(end)-BW(2);



% % bandwidth=0;X=0;
% % for Rx=60:10:990
% %     X=[X,Rx];
% %     f0=1e5;
% %     w0=2*pi*f0;
% %     f=[10:100:2e5];
% %     w=2*pi*f;
% %     V=50;
% %     Rs=50;
% %     Rl=1000;
% %     L2=((Rl-Rx)*Rx)^0.5/w0;
% %     C2=(Rl-Rx)^0.5/(Rl*Rx^0.5*w0);
% %     L1=((Rx-Rs)*Rs)^0.5/w0;
% %     C1=(Rx-Rs)^0.5/(Rx*Rs^0.5*w0);
% %     XL1=w*L1*1j;
% %     XC1=-1j./(w*C1);
% %     XL2=w*L2*1j;
% %     XC2=-1j./(w*C2);
% %     q=XL2+(XC2*Rl./(XC2+Rl));
% %     Z=Rs+XL1+((XC1.*q)./(q+XC1));
% %     I=V./Z;
% %     Io=I.*(XC1./(q+XC1)).*(XC2./(Rl+XC2));
% %     P=abs(Io.^2)*Rl;
% %     Phalf=max(P)/2;
% %     BW=0;
% %     for i=1:2000
% %         if P(i) > Phalf
% %             BW=[BW,f(i)];
% %         end
% %     end
% % bandwidth=[bandwidth,BW(end)-BW(2)];
% % end
% % plot(X,bandwidth);




% % % f0=1e5;
% % % w0=2*pi*f0;
% % % f=[10:100:2e5];
% % % w=2*pi*f;
% % % V=50;
% % % Rs=50;
% % % Rl=1000;
% % % Ll=Rl/(10*w0);
% % % Cl=10/(Rl*w0);
% % % L=((Rl-Rs)*Rs)^0.5/w0;
% % % C=(Rl-Rs)^0.5/(Rl*Rs^0.5*w0);
% % % XL1=w0*Ll*j;
% % % XCl=-j/(w0*Cl);
% % % Xl=L*w*j;
% % % Xc=-j./(C*w);
% % % Zl=XL1+Rl+XCl;
% % % Z=Rs+Xl+((Xc.*Zl)./(Zl+Xc));
% % % I=V./Z;
% % % Io=I.*Xc./(Zl+Xc);
% % % P=abs(Io.^2)*Rl;
% % % plot(f,P,'-r');
% % % Phalf=max(P)/2;
% % % BW=0;
% % % for i=1:2000
% % %     if P(i) > Phalf
% % %         BW=[BW,f(i)];
% % %     end
% % % end
% % % bandwidth=BW(end)-BW(2);



f0=1e5;
w0=2*pi*f0;
f=[10:100:2e5];
w=2*pi*f;
V=50;
Rs=50;
Rl=1000;
Ll=Rl/(10*2*pi*f0);
Cl=10/(Rl*2*pi*f0);
L=((Rl-Rs)*Rs)^0.5/w0;
C=(Rl-Rs)^0.5/(Rl*Rs^0.5*w0);
XLl=w0*Ll*j;
XCl=-j/(w0*Cl);
Xl=L*w*j;
Xc=-j./(C*w);
Zl=1/(1/XLl+1/Rl+1/XCl);
Z=Rs+Xl+((Xc*Zl)./(Zl+Xc));
I=V./Z;
Io=I.*(Xc./(Zl+Xc)).*(XLl/(XLl+(XCl*Rl)/(XCl+Rl))).*(XCl/(XCl+Rl));
P=abs(Io.^2)*Rl;
plot(f,P,'-b');
Phalf=max(P)/2;
BW=0;
for i=1:2000
    if P(i) > Phalf
        BW=[BW,f(i)];
    end
end
bandwidth=BW(end)-BW(2);


% f0=1e5;
% w0=2*pi*f0;
% f=[10:100:2e5];
% w=2*pi*f;
% V=50;
% Rs=50;
% Rl=1000;
% Ll=Rl/(10*w0);
% Cl=10/(Rl*w0);
% C = 1/(w0*(sqrt(Rs*(Rl-Rs))));
% L = (Rl/w0)*(sqrt(Rs/(Rl-Rs)));
% XLl=w0*Ll*j;
% XCl=-j/(w0*Cl);
% Xl=L*w*j;
% Xc=-j./(C*w);
% Zl=XLl+Rl+XCl;
% Z=Rs+Xc+((Xl*Zl)./(Zl+Xl));
% I=V./Z;
% Io=I.*Xl./(Zl+Xl);
% P=abs(Io.^2)*Rl;
% plot(f,P,'-g');
% Phalf=max(P)/2;
% BW=0;
% for i=1:2000
%     if P(i) > Phalf
%         BW=[BW,f(i)];
%     end
% end
% bandwidth=BW(end)-BW(2);



% % f0=1e5;
% % w0=2*pi*f0;
% % f=[10:100:2e5];
% % w=2*pi*f;
% % V=50;
% % Rs=50;
% % Rl=1000;
% % Ll=Rl/(10*2*pi*f0);
% % Cl=10/(Rl*2*pi*f0);
% % C = 1/(w0*(sqrt(Rs*(Rl-Rs))));
% % L = (Rl/w0)*(sqrt(Rs/(Rl-Rs)));
% % XLl=w0*Ll*j;
% % XCl=-j/(w0*Cl);
% % Xl=L*w*j;
% % Xc=-j./(C*w);
% % Zl=1/(1/XLl+1/Rl+1/XCl);
% % Z=Rs+Xc+((Xl*Zl)./(Zl+Xl));
% % I=V./Z;
% % Io=I.*(Xl./(Zl+Xl)).*(XLl/(XLl+(XCl*Rl)/(XCl+Rl))).*(XCl/(XCl+Rl));
% % P=abs(Io.^2)*Rl;
% % plot(f,P,'-g');
% % Phalf=max(P)/2;
% % BW=0;
% % for i=1:2000
% %     if P(i) > Phalf
% %         BW=[BW,f(i)];
% %     end
% % end
% % bandwidth=BW(end)-BW(2);



% % % V=5;
% % % Rs=50;
% % % Rl=1000;
% % % w=2*pi*1e5;
% % % 
% % % C=sqrt(Rl/Rs-1)/(Rl*w);
% % % L=Rs*sqrt(Rl/Rs-1)/w;
% % % Lmin=0.9*L;
% % % Lmax=1.1*L;
% % % Cmin=0.9*C;
% % % Cmax=1.1*C;
% % % Lx=L+(Lmax-Lmin)*(rand(1000,1)-0.5);
% % % Cx=C+(Cmax-Cmin)*(rand(1000,1)-0.5);
% % % k=1;
% % % for i=1:1000
% % %     for j=1:1000
% % %         k=k+1;
% % %         P(1,k)=(Rl*V^2)/((Rs+Rl-Rl*Lx(i)*Cx(j)*w^2)^2+(w*Lx(i)+Rs*Rl*w*Cx(j))^2);
% % %     end
% % % end
% % % hist(P,500);
% % % axis([0.1,0.14,0,90000]);
