function[]=Tammplasmon_TMM(s,t,ag,lbd)

%%%% Input consist of two or three coloumned vectors
%%% Wher the first coloumn has the wavelength the next has real part of
%%% refractive index and then the Imaginary part(if applicable)
%%%%% s and t are for the two dielectrics (SiO2 and TiO2 in our case)

%%%%% ALL DATA SHOULD BE AGIANST WAVELENGTH IN NANOMETERS

c0=3*10^8 ;
il=1000;%%%%%%%%number of points along the Wavelength axis
is=400;%%%%%%%%% number of values fo thickness of the metal film
l=linspace(800.0,1100.0,il);
n1=spline(ag(:,1)*10^3,ag(:,2),l);
K=spline(ag(:,1)*10^3,ag(:,3),l);
ns_r=3.0*ones(size(l));
% ns_r=spline(s(:,1),s(:,2),l);
ns_i=zeros(size(l));
% ns_i=spline(s(:,1),s(:,3),l);
nt_r=3.7*ones(size(l));
% nt_r=spline(t(:,1),t(:,2),l);
nt_i=0;
% nt_i=spline(t(:,1),t(:,3),l);
nt=nt_r+1j*nt_i;
ns=ns_r+1j*ns_i;
na=n1+1j*K;
n0=ones(size(l));% Surrounding medium (air in this case)

d_ag=linspace(0,100*10^-9,is);   %%% Thickness range of Ag film

[L,D]=meshgrid(l,d_ag);

n_p=50;     %%%%%%%%%%%% No. of pairs of Tio2 - SiO2 layers %%%%%%%%%%%%


i0=find(lbd,l,il); %Index of the central wavelength
di=300;
dt=lbd*10^-9/(4*nt_r(i0));  %%% Thickness of SiO2 
ds=lbd*10^-9/(4*ns_r(i0));  %%% Thickness of TiO2
n_spacer=nt_r(i0);
d_spacer=40*10.^-9;%lbd*10^-9/(4*n_spacer);
% dt=lbd*10^-9;  %%% Thickness of SiO2 
% ds=0.5*lbd*10^-9/ns_r(i0);  %%% Thickness of TiO2

kt_r=ones(is,il);
kt_i=ones(is,il);
ks_r=ones(is,il);
ks_i=ones(is,il);
k_spacer=ones(is,il);

for i=1:is
    k_spacer(i,:)=(2*pi*n_spacer)./(L(i,:)*10^-9);
    kr(i,:)=(2*pi*n1)./(L(i,:)*10^-9);
    ki(i,:)=(2*pi*K)./(L(i,:)*10^-9);
    kt_r(i,:)=(2*pi*nt_r)./(L(i,:)*10^-9);
    kt_i(i,:)=(2*pi*nt_i)./(L(i,:)*10^-9);
    ks_r(i,:)=(2*pi*ns_r)./(L(i,:)*10^-9);
    ks_i(i,:)=(2*pi*ns_i)./(L(i,:)*10^-9);
end
k_t=kt_r+1j*kt_i;  %Wave vector in TiO2
k_s=ks_r+1j*ks_i;  %Wave vector in SiO2
k_ag=kr+1j*ki;     %Wave vector in Ag

b_st=nt./ns;
b_ts=ns./nt;

b_0t=nt./n0;
b_0s=ns./n0;

b_t0=n0./nt;
b_s0=n0./ns;

b_a0=n0./na;
b_0a=na./n0;

b_sa=na./ns;
b_ta=na./nt;

b_at=nt./na;
b_as=ns./na;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_st=(1-b_st)./(1+b_st);
t_st= 2./(1+b_st);
% ((1-r_st.^2)./b_st).^(0.5);
% 2./(1+b_st);

r_ts=(1-b_ts)./(1+b_ts);
t_ts=2./(1+b_ts);

r_0t=(1-b_0t)./(1+b_0t);
t_0t=2./(1+b_0t);

r_t0=(1-b_t0)./(1+b_t0);
t_t0=2./(1+b_t0);

r_0a=(1-b_0a)./(1+b_0a);
t_0a=2./(1+b_0a);

r_a0=(1-b_a0)./(1+b_a0);
t_a0=2./(1+b_a0);

r_sa=(1-b_sa)./(1+b_sa);
t_sa=2./(1+b_sa);

r_as=(1-b_as)./(1+b_as);
t_as=2./(1+b_as);

r_at=(1-b_at)./(1+b_at);
t_at=2./(1+b_at);

r_ta=(1-b_ta)./(1+b_ta);
t_ta=2./(1+b_ta);

r_0s=(1-b_0s)./(1+b_0s);
t_0s=2./(1+b_0s);

r_s0=(1-b_s0)./(1+b_s0);
t_s0=2./(1+b_s0);

%%%%%%%%%%%%%%%%%% Transfer Matrices & Propogation Matrices%%%%%%%%%%%%%%%%

P_t=zeros(2,2,is,il);  %Propogation Matrix in TiO2
P_s=zeros(2,2,is,il);  %Propogation Matrix in SiO2
P_a=zeros(2,2,is,il);  %Propogation Matrix in Ag
T_ts=zeros(2,2,is,il);
T_st=zeros(2,2,is,il);
T_0s=zeros(2,2,is,il);
T_s0=zeros(2,2,is,il);
T_0t=zeros(2,2,is,il);
T_t0=zeros(2,2,is,il);
T_0a=zeros(2,2,is,il);
T_a0=zeros(2,2,is,il);
T_as=zeros(2,2,is,il);
T_at=zeros(2,2,is,il);
T_sa=zeros(2,2,is,il);
T_ta=zeros(2,2,is,il);
T_1N=zeros(2,2,is,il);
R=zeros(is,il);
R_0s=zeros(is,il);

P_t(1,1,:,:)=exp(-1i*kt_r*dt).*exp(kt_i*dt);
P_t(1,2,:,:)=zeros(size(L));
P_t(2,1,:,:)=zeros(size(L));
P_t(2,2,:,:)=exp(1i*kt_r*dt).*exp(-kt_i*dt);

P_t1(1,1,:,:)=exp(-1i*k_spacer*d_spacer);
P_t1(1,2,:,:)=zeros(size(L));
P_t1(2,1,:,:)=zeros(size(L));
P_t1(2,2,:,:)=exp(1i*k_spacer*d_spacer);

P_s(1,1,:,:)=exp(-1i*ks_r*ds).*exp(ks_i*ds);
P_s(1,2,:,:)=zeros(size(L));
P_s(2,1,:,:)=zeros(size(L));
P_s(2,2,:,:)=exp(1i*ks_r*ds).*exp(-ks_i*ds);

P_a(1,1,:,:)=exp(-1i*(k_ag).*D);
P_a(1,2,:,:)=zeros(size(L));
P_a(2,1,:,:)=zeros(size(L));
P_a(2,2,:,:)=exp(1i*(k_ag).*D);


for i=1:is
    T_ts(1,1,i,:)=1./t_ts;
    T_ts(1,2,i,:)=-r_st./t_ts;
    T_ts(2,1,i,:)=r_ts./t_ts;
    T_ts(2,2,i,:)=1./t_ts;

    T_st(1,1,i,:)=1./t_st;
    T_st(1,2,i,:)=-r_ts./t_st;
    T_st(2,1,i,:)=r_st./t_st;
    T_st(2,2,i,:)=1./t_st;

    T_0s(1,1,i,:)=1./t_0s;
    T_0s(1,2,i,:)=-r_s0./t_0s;
    T_0s(2,1,i,:)=r_0s./t_0s;
    T_0s(2,2,i,:)=1./t_0s;
    
    T_s0(1,1,i,:)=1./t_s0;
    T_s0(1,2,i,:)=-r_0s./t_s0;
    T_s0(2,1,i,:)=r_s0./t_s0;
    T_s0(2,2,i,:)=1./t_s0;

    T_0t(1,1,i,:)=1./t_0t;
    T_0t(1,2,i,:)=-r_t0./t_0t;
    T_0t(2,1,i,:)=r_0t./t_0t;
    T_0t(2,2,i,:)=1./t_0t;
    
    T_t0(1,1,i,:)=1./t_t0;
    T_t0(1,2,i,:)=-r_0t./t_t0;
    T_t0(2,1,i,:)=r_t0./t_t0;
    T_t0(2,2,i,:)=1./t_t0;

    T_at(1,1,i,:)=1./t_at;
    T_at(1,2,i,:)=-r_ta./t_at;
    T_at(2,1,i,:)=r_at./t_at;
    T_at(2,2,i,:)=1./t_at;
    
    T_ta(1,1,i,:)=1./t_ta;
    T_ta(1,2,i,:)=-r_at./t_ta;
    T_ta(2,1,i,:)=r_ta./t_ta;
    T_ta(2,2,i,:)=1./t_ta;

    T_as(1,1,i,:)=1./t_as;
    T_as(1,2,i,:)=-r_sa./t_as;
    T_as(2,1,i,:)=r_as./t_as;
    T_as(2,2,i,:)=1./t_as;
    
    T_sa(1,1,i,:)=1./t_sa;
    T_sa(1,2,i,:)=-r_as./t_sa;
    T_sa(2,1,i,:)=r_sa./t_sa;
    T_sa(2,2,i,:)=1./t_sa;

    T_0a(1,1,i,:)=1./t_0a;
    T_0a(1,2,i,:)=-r_a0./t_0a;
    T_0a(2,1,i,:)=r_0a./t_0a;
    T_0a(2,2,i,:)=1./t_0a;
    
    T_a0(1,1,i,:)=1./t_a0;
    T_a0(1,2,i,:)=-r_0a./t_a0;
    T_a0(2,1,i,:)=r_a0./t_a0;
    T_a0(2,2,i,:)=1./t_a0;
    
    R_0s(i,:)=r_0s;
end
    
T_pair=zeros(2,2,is,il);
T_f=zeros(2,2,is,il);
for ii= 1:il
    for jj= 1:is
        T_pair(:,:,jj,ii)=((T_ts(:,:,jj,ii)*P_s(:,:,jj,ii)*T_st(:,:,jj,ii)*P_t(:,:,jj,ii))^(n_p));
%         T_pair(:,:,jj,ii)=((P_t(:,:,jj,ii)*T_ts(:,:,jj,ii)*P_s(:,:,jj,ii)*T_st(:,:,jj,ii))^(n_p));
%         T_f(:,:,jj,ii)=T_0t(:,:,jj,ii)*T_pair(:,:,jj,ii)*P_t1(:,:,jj,ii)*T_ta(:,:,jj,ii)*P_a(:,:,jj,ii)*T_a0(:,:,jj,ii);
        T_1N(:,:,jj,ii)=P_a(:,:,jj,ii)*T_at(:,:,jj,ii)*P_t(:,:,jj,ii)*T_pair(:,:,jj,ii)*T_t0(:,:,jj,ii);
        T_f(:,:,jj,ii)=T_0a(:,:,jj,ii)*P_a(:,:,jj,ii)*T_at(:,:,jj,ii)*P_t1(:,:,jj,ii)*T_pair(:,:,jj,ii)*T_t0(:,:,jj,ii); 
    end
end
r_1n=reshape((T_1N(2,1,:,:)./T_1N(1,1,:,:)),is,il);

r=reshape(abs(T_f(2,1,:,:)./T_f(1,1,:,:)),is,il);
r_DBR=reshape((T_pair(2,1,1,:)./T_pair(1,1,1,:)),1,il);
temp=r_ta.*r_DBR.*exp(2*1i*k_spacer(1,:)*d_spacer);
temp1=r_at;
phi=angle(temp);
phi1=angle(temp1);
R=r.^2;
T=reshape(abs(1./T_f(1,1,:,:)).^2,is,il); %%%%%%%%%%%%%
% temp1=abs(T_f(2,1,:,:)./T_f(1,1,:,:)).^2;
% R=reshape(temp1,il,1);
% temp2=abs(1./T_f(1,1,:,:)).^2;
% T=reshape(temp2,is,1);
figure;
subplot(2,1,1);
y1=abs(r(150,:).^2);
% plot(l(10:il-10),y1(10:il-10))
plot(l,y1);
xlabel('Wavelength(nm)');
ylabel('Reflectivity (R)');

subplot(2,1,2); 
% y2 = angle(1i*sin(k_s*ds).*exp(1i*k_s*ds)./(1-(R_0s.^2).*exp(1i*2*k_s*ds)));%reshape(T_as(2,2,200,:),1,il);
% plot(l,y2)
% hold on
y1=reshape(phi/pi,il,1);
y2=reshape(phi1/pi,il,1);
plot(l,y1)
% hold on
% plot(l,y2)
% hold off
xlabel('Wavelength(nm)');
ylabel('Phase/\pi');
figure;
y2=reshape(phi1/pi,il,1);
plot(l,y2)
figure;
s=surf(L,D,R);
s.EdgeColor='none';
colormap 'jet';

plot_peaks(R,l);
% fileID = fopen('50-pair GaAs AlAs spacer(GaAs)300nm Sub_Air.txt','w');
% fprintf(fileID,'%6s %12s\n','l','R');
% fprintf(fileID,'%6.2f %12.8f\n',R);
% fclose(fileID);

end
function[io]=find(L,l,il)

for jj=1:il
    if(l(jj)>=L)
        io= jj-1;
        return
    end
end

end
function[]=plot_peaks(R,l)

figure;
hold on
x=[40,60,80,100,120,140,160,180,200];
for i=x
    plot(l,1-R(i,:));
end

end