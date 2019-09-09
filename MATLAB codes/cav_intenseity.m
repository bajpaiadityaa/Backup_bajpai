function[]=cav_intenseity(l0,n0,k0,L,ns0,ls)
c0=3*10^8 ;
il=1000;
is=100;

l=linspace(210.0,750.0,il);
n=spline(l0,n0,l);
K=spline(l0,k0,l);
ns=spline(ls,ns0,l);

i0=find(L,l,il);

d1=10^-9;%top silver layer
d2=155*10^-9;
% ((L*10^-9)/(ns(i0)*2))+((L*10^-9)/(ns(i0)*8));%SiO2 layer
d3=10^-6;%Bottum silver layer
kr=2*pi*n(i0)/(L*10^-9);
ki=2*pi*K(i0)/(L*10^-9);
k_ag=kr+1j*ki;
k_s=2*pi*ns(i0)/(L*(10^-9));
k0=2*pi/(L*(10^-9));
w=2*pi*c0/(L*10^-9) ;
d2=d2+pi/(2*k_s);
x=ones(is,4);
Ei=ones(410,1);
Er=zeros(410,1);

% x(:,1)=linspace(d3,0,is);
% x(:,2)=linspace(d2,0,is);
% x(:,3)=linspace(d1,0,is);
% x(:,4)=linspace(10^-6,0,is);

x(:,1)=linspace(0,d3,is);
x(:,2)=linspace(0,d2,is);
x(:,3)=linspace(0,d1,is);
x(:,4)=linspace(0,10^-6,is);
k=[k_ag,k_s,k_ag,k0];
b=[((c0*k_ag)/w).^-1,((c0*k_ag)/(ns(i0)*w)),(c0*k_ag)/(ns(i0)*w).^-1,((c0*k_ag)/w)];


for j=1:4
    a=Ei(401-100*(j-1));
    c=Er(401-100*(j-1));
    syms e1 e2 
    eq1=e1+e2==a+c;
    eq2=e1-e2==b(j)*(a-c);
    sol=solve([eq1,eq2],[e1 e2]);
    
    Ei((301-100*(j-1)):(400-100*(j-1)))=(sol.e1)*exp(-1j*k(j)*x(:,j));
    Er((301-100*(j-1)):(400-100*(j-1)))=(sol.e2)*exp(1j*k(j)*x(:,j));
   


end

d=ones(410,1);
d(1:100,1)=(10^-6)-x(:,4);
d(101:200,1)=10^-6+d1-x(:,3);
d(201:300,1)=10^-6+d1+d2-x(:,2);
d(301:400,1)=10^-6+d1+d2+d3-x(:,1);
d(401:410,1)=10^-6+d1+d2+d3+linspace(0,10^-8,10);
[x,D]=meshgrid(linspace(0,10^-8,50),d);

I=ones(410,50);
for i=1:50    
    I(:,i)=abs((Er+Ei)/Ei(1,1)).^2;
end
s=surf(x,D,I);
s.EdgeColor='none';
colormap 'jet';
figure
s=surf(x(101:300,:),D(101:300,:),I(101:300,:));
s.EdgeColor='none';
colormap 'jet';
figure
plot(d,abs((Er)./Ei).^2);
end
function[io]=find(L,l,il)

for jj=1:il
    if(l(jj)>=L)
        io= jj-1;
        return
    end
end

end