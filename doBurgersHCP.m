function [ f, A, lb,systems,bt] = doBurgersHCP(ebsd,phaseNum,poisson)
%% Preamble
%this code outputs array f (containing energies for all dislocation types)
%array A (containing XX), 
%array lb (containing lower bounds for density of each dislocation type -i.e. zero), 
%Structure systems (containing the possibleslip planes and directions and what family they belong to)

Ntypes=0;

%import the crystal symmetry from the loaded ebsd map.  iterated by calling
%function to cover all phases.  To test code type in i=1 first.
CS=ebsd(ebsd.phase==phaseNum).CS;

%get burgers vectors from crystal info, dividing by 1E4 to convert to
%microns.  magBurgPrism = length of a, magBurgPyra=sqrt(a^2+c^2)
magBurgPrism=CS.axes.y(2)/1E4;
magBurgPyra=sqrt(CS.axes.y(2)^2 + CS.axes.z(3)^2)/1E4;



%% Prism Slip System (Edge)
%b is burgers vector, n is slip plane normal

b=Miller(1,1,-2,0,CS,'uvw');
n=Miller(1,0,-1,0,CS,'hkl');

%get all equivilent vectors (b) or planes(n), and number of them (c) 
%NOTE (c) is irrelevant and will be overwritten.


[b,c] = symmetrise(b);
[n,c] = symmetrise(n,'antipodal');

%find the cases where these are perpendicular by taking the dot product.

[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));




%reload (b) and (n) with only the valid systems 
b = b(r);
n = n(c);

%save the valid slip systems and planes into a structure named 'systems'
%for prism, there should be three unique b.
systems(1).burgers=b;
systems(1).plane=n;
systems(1).name='Prism<a>';

%convert the burgers vector index from Miller-Bravais to Miller (orthagonal)
%and store in the bt double

bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

%calculate the magnitude of the burgers vector  NOTE when calculating this
%directly using ref R112, i get the same value (3a) for all.  Not so below-
%first vector is larger.

bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

%this for loop compensates for the disparity in vector lengths created from
%using the Miller system by normalizing the vectors, and overwriting the
%vector index previously in bt with a normalised version.  We want the vectors 
%in the Miller system because XX
for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPrism;
end

%calculate the line vector which is orthoganal to both the plane normal and the
%burgers vector.

for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    bt(i+Ntypes,2,:)=t.uvw/sqrt(t.u^2+t.v^2+t.w^2);
end


% for i=1:size(b,1)
%     v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
%     v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
% %     v1=v1./norm(v1);
% %     v2=v2./norm(v2);
%     t=cross(v1,v2);
%     bt(i+Ntypes,2,:)=t/norm(t);
% end


Ntypes=size(bt,1);
prismTypes=1:Ntypes;

%% Screw Dislocations <a>
%burgers vector of the screw dislocation, which is parallel to the line
%vector for screw
b=Miller(1,1,-2,0,CS,'uvw');

%find all equivalents
[b,c] = symmetrise(b);

%save the valid slip systems and planes into a structure named 'systems'
%for <a> type screw, there should be three unique b.

systems(2).burgers=b;
systems(2).plane='screw';
systems(2).name='screw<a>';

%convert this to Miller indexes
bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

%calculate the magnitude of the vector and normalise the burgers index to
%have a magnitude of <a>, the dislocation line is normalised to a value of
%1 and saved

bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

for i=1:size(b,1)
    bt(i+Ntypes,2,1:3)=bt(i+Ntypes,1,1:3)/bmag(i);
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)/bmag(i).*magBurgPrism;
end

Ntypes=size(bt,1);
screwATypes=(prismTypes(end)+1):Ntypes;

%% Basal Slip System (Edge)
%b is burgers vector, n is slip plane normal
b=Miller(1,1,-2,0,CS,'uvw');
n=Miller(0,0,0,1,CS,'hkl');

%get all equivalent vectors (b) or planes(n), and number of them (c) 
%NOTE (c) is irrelevant and will be overwritten.

[b,c] = symmetrise(b);
[n,c] = symmetrise(n,'antipodal');

%find the cases where these are perpendicular by taking the dot product (=0).
[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

%reload (b) and (n) with only the valid systems 
b = b(r);
n = n(c);

%save the valid slip systems and planes into a structure named 'systems'
%for basal slip, there should be three unique b.
systems(3).burgers=b;
systems(3).plane=n;
systems(3).name='basal<a>';

%convert the burgers vector index from Miller-Bravais to Miller (orthagonal)
%and store in the bt double
bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

%calculate the magnitude of the burgers vector  
bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

%this for loop compensates for the disparity in vector lengths created from
%using the Miller system by normalizing the vectors, and overwriting the
%vector index previously in bt with a normalised version.  We want the vectors 
%in the Miller system because of reasons.

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPrism;
end

%calculate the line vector which is orthogonal to both the plane normal and the
%burgers vector, and normalise it

for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    bt(i+Ntypes,2,:)=t.uvw/sqrt(t.u^2+t.v^2+t.w^2);
end

% for i=1:size(b,1)
%     v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
%     v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
% %     v1=v1./norm(v1);
% %     v2=v2./norm(v2);
%     t=cross(v1,v2);
%     bt(i+Ntypes,2,:)=t/norm(t);
% end

Ntypes=size(bt,1);
basalTypes=(screwATypes(end)+1):Ntypes;






%% Pyramidal Slip System (Edge)
%b is burgers vector, n is slip plane normal
b=Miller(1,1,-2,3,CS,'uvw');
n=Miller(1,0,-1,1,CS,'hkl');

%get all equivalent vectors (b) or planes(n), and number of them (c) 
%NOTE (c) is irrelevant and will be overwritten.


[b,c] = symmetrise(b);
[n,c] = symmetrise(n,'antipodal');

%find the cases where these are perpendicular by taking the dot product.
[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

%reload (b) and (n) with only the valid systems 
b = b(r);
n = n(c);

%save the valid slip systems and planes into a structure named 'systems'
%for prismatic slip, there should be SIX unique b.
systems(4).burgers=b;
systems(4).plane=n;
systems(4).name='Pyramdial<c+a>';

%convert the burgers vector index from Miller-Bravais to Miller (orthagonal)
%and store in the bt double
bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);
%calculate the magnitude of the burgers vector  
bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);
%this for loop compensates for the disparity in vector lengths created from
%using the Miller system by normalizing the vectors, and overwriting the
%vector index previously in bt with a normalised version. 

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPyra;
end

%calculate the line vector which is orthogonal to both the plane normal and the
%burgers vector, normalised to 1

for i=1:size(b,1)
    v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
    v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
%     v1=v1./norm(v1);
%     v2=v2./norm(v2);
    t=cross(v1,v2);
    bt(i+Ntypes,2,:)=t/norm(t);
end

%this is a counter for the # of systems
Ntypes=size(bt,1);
pyramidalTypes=(basalTypes(end)+1):Ntypes;


%% Pyramidal Slip System 2 (Edge)
%b is burgers vector, n is slip plane normal
b=Miller(1,1,-2,3,CS,'uvw');
n=Miller(1,1,-2,2,CS,'hkl');

%get all equivalent vectors (b) or planes(n), and number of them (c) 
%NOTE (c) is irrelevant and will be overwritten.


[b,c] = symmetrise(b);
[n,c] = symmetrise(n,'antipodal');

%find the cases where these are perpendicular by taking the dot product.
[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

%reload (b) and (n) with only the valid systems 
b = b(r);
n = n(c);

%save the valid slip systems and planes into a structure named 'systems'
%for prismatic slip, there should be SIX unique b.
systems(5).burgers=b;
systems(5).plane=n;
systems(5).name='Pyramdial2<c+a>';

%convert the burgers vector index from Miller-Bravais to Miller (orthagonal)
%and store in the bt double
bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);
%calculate the magnitude of the burgers vector  
bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);
%this for loop compensates for the disparity in vector lengths created from
%using the Miller system by normalizing the vectors, and overwriting the
%vector index previously in bt with a normalised version. 

for i=1:size(b,1)
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPyra;
end

%calculate the line vector which is orthogonal to both the plane normal and the
%burgers vector, normalised to 1


for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    bt(i+Ntypes,2,:)=t.uvw/sqrt(t.u^2+t.v^2+t.w^2);
end

% for i=1:size(b,1)
%     v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
%     v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
% %     v1=v1./norm(v1);
% %     v2=v2./norm(v2);
%     t=cross(v1,v2);
%     bt(i+Ntypes,2,:)=t/norm(t);
% end

%this is a counter for the # of systems
Ntypes=size(bt,1);
pyramidalTypes2=(pyramidalTypes(end)+1):Ntypes;





%% Screw Dislocations <c+a>
%burgers vector of the screw dislocation, which is parallel to the line
%vector for screw
b=Miller(1,1,-2,3,CS,'uvw');

%find all equivalents
[b,c] = symmetrise(b);

%save the valid slip systems and planes into a structure named 'systems'
%there should be 6 systems
systems(6).burgers=b;
systems(6).plane='screw';
systems(6).name='screw<c+a>';

%convert the burgers vector index from Miller-Bravais to Miller (orthagonal)
%and store in the bt double
bt(Ntypes+1:Ntypes+size(b,1),1,1)=round(b.u);
bt(Ntypes+1:Ntypes+size(b,1),1,2)=round(b.v);
bt(Ntypes+1:Ntypes+size(b,1),1,3)=round(b.w);

%calculate the magnitude of the burgers vector
bmag=sqrt(round(b.u).^2 + round(b.v).^2 + round(b.w).^2);

%this for loop compensates for the disparity in vector lengths created from
%using the Miller system by normalizing the vectors, and overwriting the
%vector index previously in bt with a normalised version.  

%calculate the line vector which is orthogonal to both the plane normal and the
%burgers vector.
for i=1:size(b,1)
    bt(i+Ntypes,2,1:3)=bt(i+Ntypes,1,1:3)./bmag(i);
    bt(i+Ntypes,1,1:3)=bt(i+Ntypes,1,1:3)./bmag(i).*magBurgPyra;
end

Ntypes=size(bt,1);
screwPyramidalTypes=(pyramidalTypes2(end)+1):Ntypes;
%% Calculate curvature matrix

DD(1:Ntypes) = bt(1:Ntypes,1,1);
GG(1:Ntypes) = bt(1:Ntypes,1,2);
HH(1:Ntypes) = bt(1:Ntypes,1,3);

QQ(1:Ntypes) = bt(1:Ntypes,2,2);
RR(1:Ntypes) = bt(1:Ntypes,2,3);
SS(1:Ntypes) = bt(1:Ntypes,2,1);

A = -1*[DD.*QQ, -DD.*QQ; DD.*RR, -DD.*RR; GG.*SS, -GG.*SS; GG.*RR, -GG.*RR; HH.*RR, -HH.*RR];


%poisson for zr = 0.34
%for Mg poisson=0.29
%calculate the energy of each dislocation type, normalised to that of a basal dislocation  
%If in the basal plane,this equals absolute value of a^2
%if screw=edge energy*(1-poisson's ratio)
f(prismTypes)=1;
f(screwATypes)=(1-poisson);
f(basalTypes)=1;

f(pyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2);
f(pyramidalTypes2)=(magBurgPyra^2)/(magBurgPrism^2);
f(screwPyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2)*(1-poisson);
f((Ntypes+1):(Ntypes*2))=f(1:Ntypes);
lb=zeros(Ntypes*2,1);

systems(1).indices=[prismTypes (Ntypes+prismTypes)];
systems(2).indices=[screwATypes (Ntypes+screwATypes)];
systems(3).indices=[basalTypes (Ntypes+basalTypes)];
systems(4).indices=[pyramidalTypes (Ntypes+pyramidalTypes)];
systems(5).indices=[pyramidalTypes2 (Ntypes+pyramidalTypes2)];
systems(6).indices=[screwPyramidalTypes (Ntypes+screwPyramidalTypes)];


end


