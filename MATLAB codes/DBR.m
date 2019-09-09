function[]=DBR(n1,n2,N_sub,N_0,n_p,lbd,theta,dbr,mode)

%%%% Input consist of two or three coloumned vectors
%%% Wher the first coloumn has the wavelength the next has real part of
%%% refractive index and then the Imaginary part(if applicable)
%%%%% s and t are for the two dielectrics (SiO2 and TiO2 in our case)

%%%%% ALL DATA SHOULD BE AGIANST WAVELENGTH IN NANOMETERS

c0=3*10^8 ;
il=size(dbr(:,1),1);%%%%%%%%number of points along the Wavelength axis

l=dbr(:,1);

% ns_r=spline(n1(:,1),n1(:,2),l);
% ns_i=spline(n1(:,1),n1(:,3),l);
% nt_r=spline(n2(:,1),n2(:,2),l);
% nt_i=spline(n2(:,1),n2(:,3),l);


ns_r=n1*ones(size(l));
ns_i=zeros(size(l));
nt_r=n2*ones(size(l));
nt_i=0;
nt=nt_r+1j*nt_i;
ns=ns_r+1j*ns_i;
n0=N_0*ones(size(l));% Surrounding medium incident side (air in this case)
n_sub=N_sub*ones(size(l));

i0=find(lbd,l,il); %Index of the central wavelength
dt=lbd*10^-9/(4*nt_r(i0));  %%% Thickness of SiO2 
ds=lbd*10^-9/(4*ns_r(i0));  %%% Thickness of TiO2


Cos_t=(1-((n0.*sind(theta))./nt).^2).^0.5;
Cos_s=(1-((n0.*sind(theta))./ns).^2).^0.5;
Cos_0=(1-(sind(theta)).^2).^0.5;
Cos_S=(1-((n0.*sind(theta))./n_sub).^2).^0.5;

kt=Cos_t.*(2*pi*nt)./(l*10^-9);
ks=Cos_s.*(2*pi*ns)./(l*10^-9);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Experimental Data%%%%%%%%%%%%%%%%%%%%%%%
plot(l,dbr(:,2)/100.0,'DisplayName','\bf Experiment','LineWidth',1.5);
xlabel('\bf Wavelength(nm)','FontSize',12);
ylabel('\bf Reflectivity (R)','FontSize',12);
title(['\bf SN223(TiO_2/SiO_2 5.5 bilayers on glass) incident form film side (\theta_{i} =',num2str(theta),'^0 )']);
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%Fitting Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loop=1;
l_curr=0;
while l_curr<loop
    l_curr=l_curr+1;
    if strcmp(mode,'TE')
        b_st=(nt.*Cos_t)./(ns.*Cos_s);
        b_ts=(ns.*Cos_s)./(nt.*Cos_t);

        b_0t=(nt.*Cos_t)./(n0.*Cos_0);
        b_0s=(ns.*Cos_s)./(n0.*Cos_0);

        b_t0=(n0.*Cos_0)./(nt.*Cos_t);
        b_s0=(n0.*Cos_0)./(ns.*Cos_s);

        b_sS=(n_sub.*Cos_S)./(ns.*Cos_s);
        b_tS=(n_sub.*Cos_S)./(nt.*Cos_t);

        b_Ss=(ns.*Cos_s)./(n_sub.*Cos_S);
        b_St=(nt.*Cos_t)./(n_sub.*Cos_S);

    elseif strcmp(mode,'TM')
        b_st=(ns.*Cos_t)./(nt.*Cos_s);
        b_ts=(nt.*Cos_s)./(ns.*Cos_t);

        b_0t=(n0.*Cos_t)./(nt.*Cos_0);
        b_0s=(n0.*Cos_s)./(ns.*Cos_0);

        b_t0=(nt.*Cos_0)./(n0.*Cos_t);
        b_s0=(ns.*Cos_0)./(n0.*Cos_s);

        b_sS=(ns.*Cos_S)./(n_sub.*Cos_s);
        b_tS=(nt.*Cos_S)./(n_sub.*Cos_t);

        b_Ss=(n_sub.*Cos_s)./(ns.*Cos_S);
        b_St=(n_sub.*Cos_t)./(nt.*Cos_S);
    
    elseif strcmp(mode,'TEM')
        loop=2;
        mode='TE';
        b_st=(ns.*Cos_t)./(nt.*Cos_s);
        b_ts=(nt.*Cos_s)./(ns.*Cos_t);

        b_0t=(n0.*Cos_t)./(nt.*Cos_0);
        b_0s=(n0.*Cos_s)./(ns.*Cos_0);

        b_t0=(nt.*Cos_0)./(n0.*Cos_t);
        b_s0=(ns.*Cos_0)./(n0.*Cos_s);

        b_sS=(ns.*Cos_S)./(n_sub.*Cos_s);
        b_tS=(nt.*Cos_S)./(n_sub.*Cos_t);

        b_Ss=(n_sub.*Cos_s)./(ns.*Cos_S);
        b_St=(n_sub.*Cos_t)./(nt.*Cos_S);
        
    else
        disp('Enter appropriate mode type: TE/TM/TEM');
        return
    end

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

    r_sS=(1-b_sS)./(1+b_sS);
    t_sS=2./(1+b_sS);

    r_Ss=(1-b_Ss)./(1+b_Ss);
    t_Ss=2./(1+b_Ss);

    r_St=(1-b_St)./(1+b_St);
    t_St=2./(1+b_St);

    r_tS=(1-b_tS)./(1+b_tS);
    t_tS=2./(1+b_tS);

    r_0s=(1-b_0s)./(1+b_0s);
    t_0s=2./(1+b_0s);

    r_s0=(1-b_s0)./(1+b_s0);
    t_s0=2./(1+b_s0);

    %%%%%%%%%%%%%%%%%% Transfer Matrices & Propogation Matrices%%%%%%%%%%%%%%%%

    P_t=zeros(2,2,il);  %Propogation Matrix in TiO2
    P_s=zeros(2,2,il);  %Propogation Matrix in SiO2
    T_ts=zeros(2,2,il);
    T_st=zeros(2,2,il);
    T_0s=zeros(2,2,il);
    T_s0=zeros(2,2,il);
    T_0t=zeros(2,2,il);
    T_t0=zeros(2,2,il);
    T_Sa=zeros(2,2,il);
    T_S0=zeros(2,2,il);
    T_Ss=zeros(2,2,il);
    T_St=zeros(2,2,il);
    T_sS=zeros(2,2,il);
    T_tS=zeros(2,2,il);
    T_1N=zeros(2,2,il);
    R=zeros(il);


    P_t(1,1,:)=exp(-1i*kt*dt);
    P_t(1,2,:)=zeros(size(l));
    P_t(2,1,:)=zeros(size(l));
    P_t(2,2,:)=exp(1i*kt*dt);

    P_s(1,1,:)=exp(-1i*ks*ds);
    P_s(1,2,:)=zeros(size(l));
    P_s(2,1,:)=zeros(size(l));
    P_s(2,2,:)=exp(1i*ks*ds);



    T_ts(1,1,:)=1./t_ts;
    T_ts(1,2,:)=-r_st./t_ts;
    T_ts(2,1,:)=r_ts./t_ts;
    T_ts(2,2,:)=1./t_ts;

    T_st(1,1,:)=1./t_st;
    T_st(1,2,:)=-r_ts./t_st;
    T_st(2,1,:)=r_st./t_st;
    T_st(2,2,:)=1./t_st;

    T_0s(1,1,:)=1./t_0s;
    T_0s(1,2,:)=-r_s0./t_0s;
    T_0s(2,1,:)=r_0s./t_0s;
    T_0s(2,2,:)=1./t_0s;

    T_s0(1,1,:)=1./t_s0;
    T_s0(1,2,:)=-r_0s./t_s0;
    T_s0(2,1,:)=r_s0./t_s0;
    T_s0(2,2,:)=1./t_s0;

    T_0t(1,1,:)=1./t_0t;
    T_0t(1,2,:)=-r_t0./t_0t;
    T_0t(2,1,:)=r_0t./t_0t;
    T_0t(2,2,:)=1./t_0t;

    T_t0(1,1,:)=1./t_t0;
    T_t0(1,2,:)=-r_0t./t_t0;
    T_t0(2,1,:)=r_t0./t_t0;
    T_t0(2,2,:)=1./t_t0;

    T_St(1,1,:)=1./t_St;
    T_St(1,2,:)=-r_tS./t_St;
    T_St(2,1,:)=r_St./t_St;
    T_St(2,2,:)=1./t_St;

    T_tS(1,1,:)=1./t_tS;
    T_tS(1,2,:)=-r_St./t_tS;
    T_tS(2,1,:)=r_tS./t_tS;
    T_tS(2,2,:)=1./t_tS;

    T_Ss(1,1,:)=1./t_Ss;
    T_Ss(1,2,:)=-r_sS./t_Ss;
    T_Ss(2,1,:)=r_Ss./t_Ss;
    T_Ss(2,2,:)=1./t_Ss;

    T_sS(1,1,:)=1./t_sS;
    T_sS(1,2,:)=-r_Ss./t_sS;
    T_sS(2,1,:)=r_sS./t_sS;
    T_sS(2,2,:)=1./t_sS;

    % T_0a(1,1,:)=1./t_0S;
    % T_0a(1,2,:)=-r_a0./t_0a;
    % T_0a(2,1,:)=r_0a./t_0a;
    % T_0a(2,2,:)=1./t_0a;
    % 
    % T_a0(1,1,:)=1./t_a0;
    % T_a0(1,2,:)=-r_0a./t_a0;
    % T_a0(2,1,:)=r_a0./t_a0;
    % T_a0(2,2,:)=1./t_a0;



    T_pair=zeros(2,2,il);
    T_f=zeros(2,2,il);
    for ii= 1:il
        %%% T_pair below for SN223 %%%
        T_pair(:,:,ii)=T_0t(:,:,ii)*P_t(:,:,ii)*T_ts(:,:,ii)*P_s(:,:,ii)*T_st(:,:,ii)*P_t(:,:,ii)*((T_ts(:,:,ii)*P_s(:,:,ii)*T_st(:,:,ii)*P_t(:,:,ii))^(n_p-1))*T_tS(:,:,ii);

    end

    r_DBR=reshape((T_pair(2,1,:)./T_pair(1,1,:)),il,1);
    temp=r_DBR;
    phi=angle(temp);
    R=abs(r_DBR).^2;
    % T=reshape(abs(1./T_f(1,1,:,:)).^2,is,il); %%%%%%%%%%%%%
    % temp1=abs(T_f(2,1,:,:)./T_f(1,1,:,:)).^2;
    % R=reshape(temp1,il,1);
    % temp2=abs(1./T_f(1,1,:,:)).^2;
    % T=reshape(temp2,is,1);
    y1=R;
    plot(l,y1,'DisplayName','\bf Simulation','LineWidth',2);
    hold off
    
    if(loop==2)
        hold on
    end
    % subplot(2,1,1);
    % y1=R;
    % plot(l,y1);
    % xlabel('Wavelength(nm)');
    % ylabel('Reflectivity (R)');
    % 
    % subplot(2,1,2); 
    % y1=reshape(phi/pi,il,1);
    % plot(l,y1)
    % xlabel('Wavelength(nm)');
    % ylabel('Phase/\pi');
end 
end
function[io]=find(L,l,il)

for jj=1:il
    if(l(jj)<=L)
        io= jj-1;
        return
    end
end

end
% function[]=Give_P()
% 
% end
