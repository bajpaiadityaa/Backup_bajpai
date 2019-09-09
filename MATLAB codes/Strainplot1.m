function[]=Strainplot1(i,b,l,s)
siz=max(size(s));
peak=ones(siz,1);
width=ones(siz,1);
E=ones(siz,1);
sl=max(size(l));
gaussEqn = 'a1*exp(-((x-b1)/c1)^2)+d';
% gaussEqn = 'a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2)+d';
% v=floor(3.5*sl/15):floor(sl*8.5/15);
v=1:sl;
for t=1:siz
    off=(t-1);
    strtPoint=[1,654,5,off];
%     strtPoint=[1,654,5,1,654,1,off];
    i1=((max(i(:,t)-b(:,t))^-1)*(i(:,t)-b(:,t)));%+off);
%     f1 = fit(l(v),i1(v),gaussEqn,'Start', strtPoint);
%     peak(t)=f1.b1;
%     width(t)=f1.c1;
%     E(t)=1240/f1.b1;
    plot(l(v),i1(v),'DisplayName',strcat('strain=',num2str(s(t))));
    hold on
%     plot(f1,l(v),i1(v));
%     hold on
end
xlabel('Wavelength(nm)');
ylabel('Intensity(10^3counts)');
title('PL-Spectroscopy for strained WS_2 monolayer')
grid on
% grid minor

% figure;
% plot(s,E,'-*');
% E1=polyfit(s,E,1);
% E_fit=polyval(E1,s);
% errorbar(s,E,0.0005*ones(size(s)),'-');
% hold on
% plot(s,E_fit);
% xlabel('Strain(%)');
% ylabel('Energy of the transision(nm)');
% title('PL-Spectroscopy for strained WS_2 monolayer')
% disp(['Equation is y = ' num2str(E1(1)) '*x + ' num2str(E1(2))])
% grid on
% grid minor
% 
% figure
% errorbar(s,(2.355/(2^0.5))*width,0.15*ones(size(s)),'-*');
% hold on;
% xlabel('Strain(%)');
% ylabel('Line Width (nm)');
% title('PL-Spectroscopy for strained WS_2 monolayer')
% grid on
% grid minor
% 
% figure
% errorbar(s,peak,0.15*ones(size(s)),'-*');
% hold on;
% xlabel('Strain(%)');
% ylabel('Peak wavelength (nm)');
% title('PL-Spectroscopy for strained WS_2 monolayer')
% grid on
% grid minor
end