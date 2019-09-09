function[]=interac_cav(l0,n0,k0,ns0,ls)
L=500;
il=1000;
% is=5;
l=linspace(210.0,750.0,il);
f=figure;

p=plot(l,reflec(l0,n0,k0,ns0,ls,L));

b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],'value',L, 'min',210, 'max',750);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],'String','1','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],'String','width of SiO_2','BackgroundColor',bgcolor);

b.Callback = @reflec; 

ylabel('Reflectivity');
xlabel('Wavelength(nm)');


end
function[R]=reflec(l0,n0,k0,ns0,ls,L)

c0=3*10^8 ;
il=1000;
% is=5;
l=linspace(210.0,750.0,il);
n=spline(l0,n0,l);
K=spline(l0,k0,l);
ns=spline(ls,ns0,l);
w=2*pi*c0./(l*10^-9); 

i0=find(L,l,il);
d=ones(3,1);
d(3)=10^-9;
d(2)=((L*10^-9)/(ns(i0)*2));
d(1)=10^-6;

b=ones(il,3);
k=ones(il,3);
kr=(2*pi*n)./(l*10^-9);
ki=(2*pi*K)./(l*10^-9);
k_ag=kr+1j*ki;
k_s=(2*pi*ns)./(l*(10^-9));
k(:,1)=k_ag;
k(:,2)=k_s;
k(:,3)=k_ag;

b(:,1)=((c0*k_ag)./(ns.*w)).^-1;
b(:,2)=(c0*k_ag)./(ns.*w);
b(:,3)=((c0*k_ag)./w).^-1;

A=ones(il,4);
B=ones(il,3);
A(:,1)=((1-b(:,3))./(1+b(:,3)));

for j=1:3
    B(:,j)=(b(:,j).^-1).*((1-A(:,j).*(cos(2*k(:,j).*d(j))+1j*sin(2*k(:,j).*d(j))))./(1+A(:,j).*(cos(2*k(:,j).*d(j))+1j*sin(2*k(:,j).*d(j)))));
    A(:,j+1)=(1-B(:,j))./(B(:,j)+1);
end
R=(abs(A(:,4)).^2);

end
function[io]=find(L,l,il)
for j=1:il
    if(l(j)>=L)
        io= j-1;
        return
    end
end

end