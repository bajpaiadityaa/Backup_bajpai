function[]=MfD(l0,n0,k0,L,ns0,ls)
c0=3*10^8 ;
il=10000;
is=100;

D1=linspace(10^-9,60*10^-9,is);  % Thickness of Ag film1
D2=linspace(10^-9,10^-6,is);  % Thickness of SiO2 film
[d1,d2]=meshgrid(D1,D2);

l=linspace(210.0,750.0,il);
n=spline(l0,n0,l);
K=spline(l0,k0,l);
ns=spline(ls,ns0,l);
w=2*pi*c0/(L*10^-9) ;

A=ones(is,is,3);
B=ones(is,is,2);
d=ones(is,is,2);
d(:,:,1)=d2;
d(:,:,2)=d1;

kr=(2*pi*n)./(l*10^-9);
ki=(2*pi*K)./(l*10^-9);
k_ag=kr+1j*ki;
i=find(L,l,il);
k_s=2*pi*ns(i)/(L*(10^-9));
k=[k_s,k_ag(i)];
b=[(c0*k_ag(i))/(ns(i)*w),((c0*k_ag(i))/w)^-1];
A(:,:,1)=((-ns(i)+1)/(1+ns(i)))*ones(is,is,1);
for j=1:2
    B(:,:,j)=(b(j)^-1).*(1-A(:,:,j).*(cos(2*k(j).*d(:,:,j))+1j*sin(2*k(j).*d(:,:,j))))./(1+A(:,:,j).*(cos(2*k(j).*d(:,:,j))+1j*sin(2*k(j).*d(:,:,j))));
    A(:,:,j+1)=(1-B(:,:,j))./(1+B(:,:,j));
end
R=abs(A(:,:,3)).^2;

s=surf(D1,D2,R);
s.EdgeColor='none';
grid on
grid minor
xlabel('Width of Ag film');
ylabel('Width of SiO_2 film');

title('Reflectance of Ag on SiO_2 substrate');

end

function[io]=find(L,l,il)
for j=1:il
    if(l(j)>=L)
        io= j-1;
        return
    end
end

end