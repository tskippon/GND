function [ f, A, lb,systems] = doBurgersHCP(ebsd,i,poisson)

Ntypes=0;
CS=ebsd(ebsd.phase==i).CS;

%get burgers vectors from crystal info, dividing by 1E4 to convert to
%microns
magBurgPrism=CS.axes.y(2)/1E4;
magBurgPyra=sqrt(CS.axes.y(2)^2 + CS.axes.z(3)^2)/1E4;








%Prism Slip System (Edge)
b=Miller(1,1,-2,0,CS,'uvw');
n=Miller(1,-1,0,0,CS,'hkl');


[b,c] = symmetrise(b,'antipodal');
[n,c] = symmetrise(n,'antipodal');



[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

systems(1).burgers=b;
systems(1).plane=n;
systems(1).name='Prism<a>';

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPrism;
end

for i=1:size(b,1)
    v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
    v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end


Ntypes=size(bt,1);
prismTypes=1:Ntypes;

%Screw Dislocations <a>
b=Miller(1,1,-2,0,CS,'uvw');


[b,c] = symmetrise(b,'antipodal');

systems(2).burgers=b;
systems(2).plane='screw';
systems(2).name='screw<a>';

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,2,1:3)=bt(i+Ntypes,1,1:3)/bmag(i);
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)/bmag(i).*magBurgPrism;
end

Ntypes=size(bt,1);
screwATypes=(prismTypes(end)+1):Ntypes;

%Basal Slip System (Edge)
b=Miller(1,1,-2,0,CS,'uvw');
n=Miller(0,0,0,1,CS,'hkl');



[b,c] = symmetrise(b,'antipodal');
[n,c] = symmetrise(n,'antipodal');


[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

systems(3).burgers=b;
systems(3).plane=n;
systems(3).name='basal<a>';

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPrism;
end

for i=1:size(b,1)
    v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
    v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end

Ntypes=size(bt,1);
basalTypes=(screwATypes(end)+1):Ntypes;







%Pyramidal Slip System (Edge)

b=Miller(1,1,-2,3,CS,'uvw');
n=Miller(1,0,-1,1,CS,'hkl');


[b,c] = symmetrise(b,'antipodal');
[n,c] = symmetrise(n,'antipodal');



[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

systems(4).burgers=b;
systems(4).plane=n;
systems(4).name='Pyramdial<c+a>';


bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPyra;
end

for i=1:size(b,1)
    v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
    v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end

Ntypes=size(bt,1);
pyramidalTypes=(basalTypes(end)+1):Ntypes;


%Screw Dislocations <c+a>
b=Miller(1,1,-2,3,CS,'uvw');


[b,c] = symmetrise(b,'antipodal');

systems(5).burgers=b;
systems(5).plane='screw';
systems(5).name='screw<c+a>';

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,2,1:3)=bt(i+Ntypes,1,1:3)./bmag(i);
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPyra;
end

Ntypes=size(bt,1);
screwPyramidalTypes=(pyramidalTypes(end)+1):Ntypes;

DD(1:Ntypes) = bt(1:Ntypes,1,1);
GG(1:Ntypes) = bt(1:Ntypes,1,2);
HH(1:Ntypes) = bt(1:Ntypes,1,3);

QQ(1:Ntypes) = bt(1:Ntypes,2,2);
RR(1:Ntypes) = bt(1:Ntypes,2,3);
SS(1:Ntypes) = bt(1:Ntypes,2,1);

A = -1*[DD.*QQ, -DD.*QQ; DD.*RR, -DD.*RR; GG.*SS, -GG.*SS; GG.*RR, -GG.*RR; HH.*RR, -HH.*RR];


%poisson for zr = 0.34
f(prismTypes)=1;
f(screwATypes)=(1-poisson);
f(basalTypes)=1;
f(pyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2);
f(screwPyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2)*(1-poisson);
f((Ntypes+1):(Ntypes*2))=f(1:Ntypes);
lb=zeros(Ntypes*2,1);

systems(1).indices=[prismTypes (Ntypes+prismTypes)];
systems(2).indices=[screwATypes (Ntypes+screwATypes)];
systems(3).indices=[basalTypes (Ntypes+basalTypes)];
systems(4).indices=[pyramidalTypes (Ntypes+pyramidalTypes)];
systems(5).indices=[screwPyramidalTypes (Ntypes+screwPyramidalTypes)];

end


