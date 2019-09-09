function[]=Metal_film(n_m,l_m,k_m,L)
c0=3*10^8 ;
il=1000;%%%%%%%%number of points along the Wavelength axis
is=100;%%%%%%%%% number of values fo thickness of the metal film
r=ones(is,1);
l=linspace(210,750,il);
s=linspace(10^-9,70*10^-9,is);
n=spline(l_m,n_m,l);
K=spline(l_m,k_m,l);
w=2*pi*c0./(l*10^-9) ;
%%%%%%% Real and imaginery part of the K vector of the metal film
kr=(2*pi.*n)./(l*10^-9);
ki=(2*pi.*K)./(l*10^-9);
k=kr+1j*ki;
i0=find(L,l,il);
a=((n(i0)-1)+1j*K(i0))/((n(i0)+1)+1j*K(i0));

% h=1-(((n(i0)-1).^2-(k(i0)^2)+2*(n(i0)-1)*K(i0)*1j)./((n(i0)+1)^2-(k(i0)^2)+2*(n(i0)+1)*K(i0)*1j)).*(exp(-2*ki(i0)*s).*(cos(2*kr(i0)*s)+1j*sin(2*kr(i0)*s)));
o=(abs(a).^2);
r=o.*abs((1-exp(-2*ki(i0)*s).*(cos(2*kr(i0)*s)+1j*sin(2*kr(i0)*s)))).^2;

plot(s,r);
xlabel('Thickness');
ylabel('Reflectivity');

title('Reflectivity of Silver sheet');
grid on
grid minor
end

function[io]=find(L,l,il)
for j=1:il
    if(l(j)>=L)
        io= j-1;
        return
    end
end

end