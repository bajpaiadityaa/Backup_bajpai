function[]=MD(l_m,n_m,k_m,L,n_d,l_d)
%%%%%%%%%%% this is a simulation of the reflectivity spectra of a metal
%%%%%%%%%%% film deposited over an infinite dieletric surface

c0=3*10^8 ;
il=10000;%%%%%% number of points on the Wavelength axis
is=100;%%%%%%% number of values of metal film thickness


s= linspace(10^-9,200*10^-9,is) ;%%%%%range of values of metal film thickness
l=linspace(210.0,750.0,il);
n=spline(l_,n_m,l);
K=spline(l_m,k_m,l);
ns=spline(l_d,n_d,l);
w=2*pi*c0./(l*10^-9) ;

%%%%%%% real and imaginary part of the K vector for metal
kr=(2*pi.*n)./(l*10^-9);
ki=(2*pi.*K)./(l*10^-9);
k=kr+1j*ki;

i=find(L,l,il);

a0=((n(i)-1)+1j*K(i))/((n(i)+1)+1j*K(i));
a1=(((n(i)/ns(i))-1)+1j*(K(i)/ns(i)))/(((n(i)/ns(i))+1)+1j*(K(i)/ns(i)));
t=exp(-2*ki(i)*s).*(cos(2*kr(i)*s)+1j*sin(2*kr(i)*s));
r1=(a1-a0*t)./(1-a0*a1*t);
r=abs(r1).^2;
plot(s*10^9,r);
xlabel('Thickness of Ag film (nm)');
ylabel('Reflectivity');
title(strcat('Reflectivity vs film thickness for a silver film on a SiO2 substarte'))
end


function[io]=find(L,l,il)
for j=1:il
    if(l(j)>=L)
        io= j-1;
        return
    end
end

end