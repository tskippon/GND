function [ f, A, lb,systems] = doBurgersHCP(ebsd,phaseNum,poisson)
%% Preamble
%this code outputs array f (containing energies for all dislocation types)
%array A (containing Burgers and Line vector information), 
%array lb (containing lower bounds for density of each dislocation type -i.e. zero), 
%Structure systems (containing the possibleslip planes and directions and what family they belong to)

Ntypes=0;

%import the crystal symmetry from the loaded ebsd map.  iterated by calling
%function to cover all phases.  To test code type in i=1 first.
CS=ebsd(ebsd.phase==phaseNum).CS;



%%Get units of ebsd coordinates and set burgers vector unit conversion
%%appropriately. (units of Miller type objects are in Angstroms by default,
%%so need to convert from that).

if strcmp(ebsd.scanUnit,'nm')
    unitConversion=1e-1;
elseif strcmp(ebsd.scanUnit,'um')
    unitConversion=1e-4;
elseif strcmp(ebsd.scanUnit,'mm')
    unitConversion=1e-7;
elseif strcmp(ebsd.scanUnit,'m')
    unitConversion=1e-10;
else
    disp('Warning! Units of EBSD scan coordinates not recognized! Assuming scan is in microns.')
    unitConversion=1e-4;
end



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

burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;



%calculate the line vector which is orthoganal to both the plane normal and the
%burgers vector.

for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    line(i+Ntypes)=t.normalize;
end


% for i=1:size(b,1)
%     v1=[round(b(i).u) round(b(i).v) round(b(i).w)];
%     v2=[round(n(i).h) round(n(i).k) round(n(i).l)];
% %     v1=v1./norm(v1);
% %     v2=v2./norm(v2);
%     t=cross(v1,v2);
%     bt(i+Ntypes,2,:)=t/norm(t);
% end


Ntypes=size(burgers,2);
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
burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;
line(Ntypes+1:Ntypes+size(b,1))=b.normalize;



Ntypes=size(burgers,2);
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
burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;


%calculate the line vector which is orthogonal to both the plane normal and the
%burgers vector, and normalise it

for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    line(i+Ntypes)=t.normalize;
end



Ntypes=size(burgers,2);
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
systems(4).name='Pyramidal<c+a>';

%convert the burgers vector index from Miller-Bravais to Miller (orthagonal)
%and store in the bt double
burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;


%calculate the line vector which is orthogonal to both the plane normal and the
%burgers vector, normalised to 1

for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    line(i+Ntypes)=t.normalize;
end

%this is a counter for the # of systems
Ntypes=size(burgers,2);
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
systems(5).name='Pyramidal2<c+a>';

%convert the burgers vector index from Miller-Bravais to Miller (orthagonal)
%and store in the bt double
burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;

%calculate the line vector which is orthogonal to both the plane normal and the
%burgers vector, normalised to 1


for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    line(i+Ntypes)=t.normalize;
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
Ntypes=size(burgers,2);
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
burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;
line(Ntypes+1:Ntypes+size(b,1))=b.normalize;


Ntypes=size(burgers,2);
screwPyramidalTypes=(pyramidalTypes2(end)+1):Ntypes;
%% Calculate curvature matrix
%assign to temporary variables
bs=burgers;
ls=line;
%Calculate A matrix
A = -1*[bs.x.*ls.x-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.x.*bs.y; ls.y.*bs.x; ls.y.*bs.y-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.z.*bs.x; ls.z.*bs.y];





%poisson for zr = 0.34
%for Mg poisson=0.29
%calculate the energy of each dislocation type, normalised to that of a basal dislocation  
%If in the basal plane,this equals absolute value of a^2
%if screw=edge energy*(1-poisson's ratio)
f(prismTypes)=1;
f(screwATypes)=(1-poisson);
f(basalTypes)=1;

magBurgPrism=norm(bs(prismTypes(1)));
magBurgPyra=norm(bs(pyramidalTypes(1)));

f(pyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2);
f(pyramidalTypes2)=(magBurgPyra^2)/(magBurgPrism^2);
f(screwPyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2)*(1-poisson);
lb=zeros(Ntypes,1);

systems(1).indices=prismTypes;
systems(2).indices=screwATypes;
systems(3).indices=basalTypes;
systems(4).indices=pyramidalTypes;
systems(5).indices=pyramidalTypes2;
systems(6).indices=screwPyramidalTypes;


end


