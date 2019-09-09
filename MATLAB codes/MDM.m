function[]=MDM(l_m,n_m,k_m,L,n_d,l_d)
%%%% this is a simulation of a finite dielectric sandviched between two
%%%% metal films (cavity)L  being the wavelength for which the reflectivity
%%%% will be calculated

c0=3*10^8 ;
il=100;%%%%%% number of points on the Wavelength axis
is=100;%%%%%%% number of values of metal film thickness

D1=linspace(10^-9,60*10^-9,is);  % Thickness of Ag film1
D2=linspace(10^-9,10^-6,is);  % Thickness of SiO2 film
D3=linspace(10^-9,60*10^-9,is);  % Thickness of Ag film2
[d1,d3,d2]=meshgrid(D1,D3,D2);

l=linspace(210.0,750.0,il);
n=spline(l_m,n_m,l); 

K=spline(l_m,k_m,l);
ns=spline(l_d,n_d,l);
w=2*pi*c0/(L*10^-9) ;

A=ones(is,is,is,4);
B=ones(is,is,is,3);
d=ones(is,is,is,3);
d(:,:,:,1)=d3;
d(:,:,:,2)=d2;
d(:,:,:,3)=d1;

kr=(2*pi*n)./(l*10^-9);
ki=(2*pi*K)./(l*10^-9);
k_ag=kr+1j*ki;
i=find(L,l,il);
k_s=2*pi*ns(i)/(L*(10^-9));
k=[k_ag(i),k_s,k_ag(i)];
b=[((c0*k_ag(i))/(ns(i)*w))^-1,(c0*k_ag(i))/(ns(i)*w),((c0*k_ag(i))/w)^-1];
A(:,:,:,1)=((1-b(3))/(1+b(3)))*ones(is,is,is,1);

for j=1:3
    B(:,:,:,j)=(b(j)^-1).*((1-A(:,:,:,j).*(cos(2*k(j).*d(:,:,:,j))+1j*sin(2*k(j).*d(:,:,:,j))))./(1+A(:,:,:,j).*(cos(2*k(j).*d(:,:,:,j))+1j*sin(2*k(j).*d(:,:,:,j)))));
    A(:,:,:,j+1)=(1-B(:,:,:,j))./(B(:,:,:,j)+1);
end
R=abs(A(:,:,:,4)).^2;
h=slice(d1,d3,d2,R,[],linspace(10^-9,60*10^-9,3),linspace(10^-9,10^-6,3));
set(h,'edgecolor','none')
grid on
grid minor
xlabel('Width of Ag film (top)');
ylabel('Width of Ag film (bottum)');
zlabel('Width of SiO2 film');
title('Reflectance dependance of SiO_2-Ag cavity');

figure
h=slice(d1,d3,d2,R,[],linspace(10^-9,60*10^-9,5),[]);
set(h,'edgecolor','none')
contourslice(d1,d3,d2,R,[],linspace(10^-9,60*10^-9,5),[]);
grid on
grid minor
xlabel('Width of Ag film (top)');
ylabel('Width of Ag film (bottum)');
zlabel('Width of SiO2 film');
title('Reflectance dependance of SiO_2-Ag cavity');

end

function[io]=find(L,l,il)
for j=1:il
    if(l(j)>=L)
        io= j-1;
        return
    end
end

end