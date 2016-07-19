function [ f, A, lb, systems] = doBurgersFCC(ebsd,i,poisson)
Ntypes=0;
CS=ebsd(ebsd.phase==1).CS;
magBurg=1/sqrt(2)*CS.axes.x(1)/1E4; %burgers vector magnitude (micrometers) NEED TO SET FOR FCC MATERIAL



%**************************
%First Slip System (Edge)
%**************************
b=Miller(1,1,0,CS,'hkl');
n=Miller(1,1,1,CS,'uvw');


[b,c] = symmetrise(b,'antipodal');
[n,c] = symmetrise(n,'antipodal');



%m = symmetrise(b);
%n = symmetrise(n);

[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

systems(1).burgers=b;
systems(1).plane=n;

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.h);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.k);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.l);

bmag=sqrt((round(b.h)).^2 + (round(b.k)).^2 + round(b.l).^2);

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurg;
end

for i=1:size(b,1)
    v1=[round(b(i).h) round(b(i).k) round(b(i).l)];
    v2=[round(n(i).U) round(n(i).V) round(n(i).W)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end

Ntypes=size(bt,1);
edge=1:Ntypes;


%**************************
%Screw Dislocations 
%**************************
b=Miller(1,1,0,CS,'hkl');
[b,c] = symmetrise(b,'antipodal');

systems(2).burgers=b;
systems(2).plane='screw';

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.h);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.k);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.l);

bmag=sqrt((round(b.h)).^2 + (round(b.k)).^2 + round(b.l).^2);
for i=1:size(b,1)
    bt(i+Ntypes,2,1:3)=bt(i+Ntypes,1,1:3)/bmag(i);
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)/bmag(i).*magBurg;
end

Ntypes=size(bt,1);
screw=(edge(end)+1):Ntypes;


DD(1:Ntypes) = bt(1:Ntypes,1,1);
GG(1:Ntypes) = bt(1:Ntypes,1,2);
HH(1:Ntypes) = bt(1:Ntypes,1,3);

QQ(1:Ntypes) = bt(1:Ntypes,2,2);
RR(1:Ntypes) = bt(1:Ntypes,2,3);
SS(1:Ntypes) = bt(1:Ntypes,2,1);

A = -1*[DD.*QQ, -DD.*QQ; DD.*RR, -DD.*RR; GG.*SS, -GG.*SS; GG.*RR, -GG.*RR; HH.*RR, -HH.*RR];

f(edge)=1;
f(screw)=1-poisson;  
f((Ntypes+1):(2*Ntypes))=f(1:Ntypes);
lb=zeros(Ntypes*2,1);

systems(1).indices=[edge (Ntypes+edge)];
systems(2).indices=[screw (Ntypes+screw)];

end