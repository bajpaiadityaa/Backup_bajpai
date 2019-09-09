function[]=tamplasmon_cav(s,t,ag,lbd)   

%%%% Input consist of two or three coloumned vectors
%%% Wher the first coloumn has the wavelength the next has real part of
%%% refractive index and then the Imaginary part(if applicable)
%%%%%s and t are for the two dielectrics (SiO2 and TiO2 in our case)
c0=3*10^8 ;
il=1000;%%%%%%%%number of points along the Wavelength axis
is=400;%%%%%%%%% number of values fo thickness of the metal film
l=linspace(800.0,1100.0,il);
n1=spline(ag(:,1)*10^3,ag(:,2),l);
K=spline(ag(:,1)*10^3,ag(:,3),l);
ns_r=spline(s(:,1),s(:,2),l);%nm data
ns_i=spline(s(:,1),s(:,3),l);
nt_r=spline(t(:,1),t(:,2),l);
nt_i=spline(t(:,1),t(:,3),l);
% ns_r=spline(s(:,1)*10^3,s(:,2),l);%micro meter data
% ns_i=0;
% ns_i=spline(s(:,1)*10^3,s(:,3),l);
% nt_r=spline(t(:,1)*10^3,t(:,2),l);
% nt_i=spline(t(:,1)*10^3,t(:,3),l);
nt=nt_r+1j*nt_i;
ns=ns_r+1j*ns_i;

i0=find(lbd,l,il);

dt=lbd*10^-9/(4*nt_r(i0));  %%% Thickness of SiO2 
ds=lbd*10^-9/(4*ns_r(i0));  %%% Thickness of TiO2
n=10;  %%% Pairs of TiO2 and SiO2 in the DBR

d_ag=linspace(10^-9,100*10^-9,is);   %%% Thickness range of Ag film

[L,D]=meshgrid(l,d_ag);

b=ones(is,il,2);
k=ones(is,il,2);
kr=ones(is,il);
ki=ones(is,il);
k_s=ones(is,il);
w=2*pi*c0./(L*10^-9);
kt_r=ones(is,il);
kt_i=ones(is,il);
ks_r=ones(is,il);
ks_i=ones(is,il);
for i=1:is
    kr(i,1:il)=(2*pi*n1)./(L(i,1:il)*10^-9);
    ki(i,1:il)=(2*pi*K)./(L(i,1:il)*10^-9);
    k_s(i,1:il)=(2*pi*ns)./(L(i,1:il)*(10^-9));
    kt_r(i,1:il)=(2*pi*nt_r)./(L(i,1:il)*10^-9);
    kt_i(i,1:il)=(2*pi*nt_i)./(L(i,1:il)*10^-9);
    ks_r(i,1:il)=(2*pi*ns_r)./(L(i,1:il)*10^-9);
    ks_i(i,1:il)=(2*pi*ns_i)./(L(i,1:il)*10^-9);
end
k_t=kt_r+1j*kt_i;
k_s=ks_r+1j*ks_i;
k_ag=kr+1j*ki;
k(:,:,1)=k_t;
k(:,:,2)=k_s;

A=ones(is,il,2*n+1);
B=ones(is,il);
b_ag=ones(is,il);
for i =1:is
    b(i,:,1)=nt./ns;
    b(i,:,2)=ns./nt;
    A(i,:,1)=((1-nt.^-1)./(1+nt.^-1));
    b_ag(i,:)=ns.*w(i,:)./(k_ag(i,:)*c0);
end

d=ones(is,il,2);
d(:,:,1)=dt*ones(is,il);
d(:,:,2)=ds*ones(is,il);
for j=1:2*n
    
   if(j~=2*n)
        B(:,:,j)=(b(:,:,mod(j,2)+1)).*((1-A(:,:,j).*(exp(1j*2*k(:,:,mod(j,2)+1).*d(:,:,mod(j,2)+1))))./(1+A(:,:,j).*(exp(1j*2*k(:,:,mod(j,2)+1).*d(:,:,mod(j,2)+1)))));
        A(:,:,j+1)=(1-B(:,:,j))./(B(:,:,j)+1);
   
   else

        B(:,:,j)=b_ag.*((1-A(:,:,j).*(exp(1j*2*k(:,:,mod(j,2)+1).*d(:,:,mod(j,2)+1))))./(1+A(:,:,j).*(exp(1j*2*k(:,:,mod(j,2)+1).*d(:,:,mod(j,2)+1)))));
        A(:,:,j+1)=(1-B(:,:,j))./(B(:,:,j)+1);
   end
end
A1=A(:,:,2*n+1);
B1=(k_ag*c0./w).*((1-A1.*(exp(1j*2*k_ag.*D)))./(1+A1.*(exp(1j*2*k_ag.*D))));
r=(1-B1)./(B1+1);
Phi_sing=angle((A1.*(exp(1j*2*k_ag.*D))));
R=abs(r).^2;
% li=il/2;
% lf=il*5/8;
% j0=is/5;
% jf=400;
% 
% r_peak=ones(jf-j0,1);
% l_peak=ones(jf-j0,1);
% Q_fac=ones(jf-j0,1);
% ref=ones(jf-j0,1);
% 
% for i= j0+1:jf  
%     inv_R=ones(1,floor(lf)-floor(li)+1)-R(i,floor(li):floor(lf));
%     l1=l(floor(li):floor(lf));
%     if(mod(i,20)==0)
%         plot(l1,inv_R,'DisplayName',num2str(d_ag(i)));
%         hold on
%     end
%     
%     
%      [max2,maxidx]=findpeaks(inv_R);
%     
%     
%     r_peak(i-j0)=inv_R(max(maxidx));   
%     l_peak(i-j0)=l1(max(maxidx));
%     [df(i-j0),Q_fac(i-j0)]=fwhm(inv_R,l1,i);
%   [df(i-j0),Q_fac(i-j0)]=fwhm(inv_R,l1,i);%%%% use either of these lines
%   depending on the function fwhm that you are using below
%     Q_fac(i-j0)=abs(l1(i1)-l1(i2));%l1(floor((i1+i2)/2))/
%     ref(i-j0)=R(max(maxidx));
% end
% plot(l_peak,r_peak*100);
% ylabel('Reflectivity');
% xlabel('Wavelength(nm)');
% figure
figure
s=surf(D,L,R);
s.EdgeColor='none';
colormap 'jet';
xlabel('thickness of Ag');
ylabel('Wavelength(nm)');
title('Reflectivity Spectra');
figure;
subplot(2,1,1);
y1=reshape(R(200,:),il,1);
plot(l,y1)
xlabel('Wavelength(nm)');
ylabel('Reflectivity (R)');

subplot(2,1,2); 
y2 = reshape(Phi_sing(200,:),il,1);
plot(l,y2)
xlabel('Wavelength(nm)');
ylabel('Phase');
set(gca,'XTick',-pi:pi/2:pi) 
set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% reflec=R(floor(0.4*is),:);
% figure;

% plot(l,reflec)
% xlabel('Wavelength');
% ylabel('Reflectivity');
% title('Reflectivity Spectra');

% figure
% plot(l_peak,Q_fac);
% ylabel('Q-factor');
% xlabel('Wavelength(nm)');
% figure
% plot(d_ag(j0+1:jf),Q_fac);
% ylabel('Q-factor');
% xlabel('Thickness of top Ag film');
% figure
% plot(Q_fac,ref*100,'*');
% xlabel('Q-factor');
% ylabel('Reflectivity');
end
function[io]=find(L,l,il)

for jj=1:il
    if(l(jj)>=L)
        io= jj-1;
        return
    end
end

end

%%%% the below commented functin calculates Qfac and delta-lambda by fiting
%%%% a gaussian

% function[df,q]=fwhm(R,l,ii)
% gaussEqn = 'a1*exp(-((x-b1)/c1)^2)+d';
% % if(ii>160)
% %     strtPoint=[1,590,5,0];
% % else
% %     strtPoint=[1,618,5,0];
% % end
% strtPoint=[1,647,5,0];
% f1 = fit(l',R',gaussEqn,'Start', strtPoint);
% peak=f1.b1;
% df=2.355*f1.c1/(2^0.5);
% 
% q=peak/df;
% if(mod(ii,20)==0)
%     plot(f1,l,R);
%     hold on
% end
% end


%%%%%%%%%%% This function calculates the fwhm using its direct definition
function[i1,i2]=fwhm(R,max,maxIdx,i_lim)
flag1=1;
flag2=1;
i1=-1;
i2=i1;
for i=1:i_lim-maxIdx-1
   if(flag1)
      if(R(i_lim)<=max/2)
         i1=maxIdx+i;
         flag1=0;
      end
   end
end
for i=0:maxIdx-1
   if(flag2)
      if(R(maxIdx-i)<=max/2)
           i2=maxIdx-i;
           flag2=0;
      end
   end
end
if(i1==-1)
    i1=1;
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
if(i2==-1)
    i2=i_lim;
end


end