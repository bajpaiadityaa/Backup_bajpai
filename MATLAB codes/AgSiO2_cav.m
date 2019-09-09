function[]=AgSiO2_cav(l_ag,n_ag,k_ag,ns0,ls,L)
c0=3*10^8 ;
il=1000;%%%%%% number of divisions on the wavelength axis
is=100;%%%%%%% number of values of Ag film width
l=linspace(200.0,1000.0,il);%%%%%%%%%% Wavelength(nm) choose the limits based on the data you have
n=spline(l_ag,n_ag,l);
K=spline(l_ag,k_ag,l);
ns=spline(ls,ns0,l);
w=2*pi*c0./(l*10^-9); %%%%angular frequency

Q_fac=ones(is,1);
ref=ones(is,1);
i0=find(L,l,il);
d=ones(3,1);
d1=linspace(10^-9,80*10^-9,is);%%%%%%%%%%Range of values of Ag film width 
d(2)= ((L*10^-9)/(ns(i0)*4))+((L*10^-9)/(ns(i0)*7.2));%%%%%%%Thickness of SiO2 film
d(1)=1;%%%%%%Thickness of the bottum Ag film

b=ones(il,3);
k=ones(il,3);
kr=(2*pi*n)./(l*10^-9);%%%%%%%%real component of the wave vector for Ag
ki=(2*pi*K)./(l*10^-9);%%%%%%%%complex component of the wave vector for Ag
k_ag=kr+1j*ki;

k_s=(2*pi*ns)./(l*(10^-9));
k(:,1)=k_ag;
k(:,2)=k_s;
k(:,3)=k_ag;

%%%% Soem constants that will be required in the calculation
b(:,1)=((c0*k_ag)./(ns.*w)).^-1;
b(:,2)=(c0*k_ag)./(ns.*w);
b(:,3)=((c0*k_ag)./w).^-1;

r_peak=ones(is,1);%%%%%%these will be used to store the lambda where the peak occours and the corresponding reflectivity
l_peak=ones(is,1);

A=ones(il,4);
B=ones(il,3);
A(:,1)=((1-b(:,3))./(1+b(:,3)));%%%% Ratio of incident and reflected E-field at the bottum most interface
for i=1:is
    d(3)=d1(i);
    for j=1:3
        B(:,j)=(b(:,j).^-1).*((1-A(:,j).*(exp(1j*2*k(:,j).*d(j))))./(1+A(:,j).*(exp(1j*2*k(:,j).*d(j)))));
        A(:,j+1)=(1-B(:,j))./(B(:,j)+1);
    end
    %%%%%%%% calculating Q-factors
    R=abs(A(:,4)).^2;
    [max_1,maxidx]=findpeaks(R);
    
    inv_R=R(max(maxidx))*ones(il,1)-R;
    [max2,maxidx]=findpeaks(inv_R);
    r_peak(i)=R(max(maxidx));   
    l_peak(i)=l(max(maxidx));
    [i1,i2]=fwhm(inv_R,inv_R(max(maxidx)),max(maxidx),il);
  
    Q_fac(i)=l(max(maxidx))/abs(l(i1)-l(i2));
    ref(i)=R(max(maxidx));
    
    if(i==floor(4.166*is/8))
        plot(l,R,'DisplayName',num2str(d(3)));
        hold on
    end
end

ylabel('Reflectivity');
xlabel('Wavelength(nm)');
plot(l_peak,r_peak*100);
ylabel('Reflectivity');
xlabel('Wavelength(nm)');
figure
plot(l_peak,Q_fac);
ylabel('Q-factor');
xlabel('Wavelength(nm)');
figure
plot(d1,Q_fac);
ylabel('Q-factor');
xlabel('Thickness of top Ag film');
figure
plot(Q_fac,ref*100,'*');
xlabel('Q-factor');
ylabel('Reflectivity');
end

function[io]=find(L,l,il)
for j=1:il
    if(l(j)>=L)
        io= j-1;
        return
    end
end

end


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