function [ f, A, lb, systems] = doBurgersBCC(ebsd,i,poisson)
Ntypes=0;
CS=ebsd(ebsd.phase==i).CS;
magBurg=(sqrt(3)/2)*CS.axes.x(1)/1E4; %burgers vector magnitude (micrometers)



%**************************
%First Slip System (Edge)
%**************************
b=Miller(1,1,1,CS,'uvw');
n=Miller(1,1,0,CS,'hkl');

[b,c] = symmetrise(b,'antipodal');
[n,c] = symmetrise(n,'antipodal');


%m = symmetrise(b);
%n = symmetrise(n);

[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

systems(1).burgers=b;
systems(1).plane=n;

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt((round(b.u)).^2 + (round(b.v)).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurg;
end

for i=1:size(b,1)
    v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
    v2=[round(n(i).k) round(n(i).k) round(n(i).l)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end

Ntypes=size(bt,1);
first=1:Ntypes;

%**************************
%Second Slip System (Edge)
%**************************
b=Miller(1,1,1,CS,'uvw');
n=Miller(1,1,2,CS,'hkl');


[b,c] = symmetrise(b,'antipodal');
[n,c] = symmetrise(n,'antipodal');

systems(2).burgers=b;
systems(2).plane=n;

%m = symmetrise(b);
%n = symmetrise(n);

[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt((round(b.u)).^2 + (round(b.v)).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurg;
end

for i=1:size(b,1)
    v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
    v2=[round(n(i).k) round(n(i).k) round(n(i).l)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end

Ntypes=size(bt,1);
secnd=(first(end)+1):Ntypes;

%**************************
%Third Slip System (Edge)
%**************************
b=Miller(1,1,1,CS,'uvw');
n=Miller(1,2,3,CS,'hkl');


[b,c] = symmetrise(b,'antipodal');
[n,c] = symmetrise(n,'antipodal');

systems(3).burgers=b;
systems(3).plane=n;

%m = symmetrise(b);
%n = symmetrise(n);

[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt((round(b.u)).^2 + (round(b.v)).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurg;
end

for i=1:size(b,1)
    v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
    v2=[round(n(i).k) round(n(i).k) round(n(i).l)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end

Ntypes=size(bt,1);
third=(secnd(end)+1):Ntypes;
%**************************
%Screw Dislocations 
%**************************
b=Miller(1,1,1,CS,'uvw');
[b,c] = symmetrise(b,'antipodal');

systems(4).burgers=b;
systems(4).plane='screw';

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt((round(b.u)).^2 + (round(b.v)).^2 + round(b.w).^2);
for i=1:size(b,1)
    bt(i+Ntypes,2,1:3)=bt(i+Ntypes,1,1:3)/bmag(i);
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)/bmag(i).*magBurg;
end

Ntypes=size(bt,1);
screw=(third(end)+1):Ntypes;





DD(1:Ntypes) = bt(1:Ntypes,1,1);
GG(1:Ntypes) = bt(1:Ntypes,1,2);
HH(1:Ntypes) = bt(1:Ntypes,1,3);

QQ(1:Ntypes) = bt(1:Ntypes,2,2);
RR(1:Ntypes) = bt(1:Ntypes,2,3);
SS(1:Ntypes) = bt(1:Ntypes,2,1);

A = -1*[DD.*QQ, -DD.*QQ; DD.*RR, -DD.*RR; GG.*SS, -GG.*SS; GG.*RR, -GG.*RR; HH.*RR, -HH.*RR];

f([first secnd third])=1;
f(screw)=1-poisson;  %Poisson ratio for beta Zr =0.45
f((Ntypes+1):(Ntypes*2))=f(1:Ntypes);
lb=zeros(Ntypes*2,1);


systems(1).indices=[first (Ntypes+first)];
systems(2).indices=[secnd (Ntypes+secnd)];
systems(3).indices=[third (Ntypes+third)];
systems(4).indices=[screw (Ntypes+screw)];
end