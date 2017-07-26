function [ f, A, lb, systems] = doBurgersFCC(ebsd,phaseNum,poisson)
Ntypes=0;
CS=ebsd(ebsd.phase==phaseNum).CS;


%%Get units of ebsd coordinates and set burgers vector unit conversion
%%appropriately.  (units of Miller type objects are in Angstroms by default,
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


%**************************
%First Slip System (Edge)
%**************************
b=Miller(1,1,0,CS,'uvw');
n=Miller(1,1,1,CS,'hkl');


[b,c] = symmetrise(b);
[n,c] = symmetrise(n,'antipodal');



[r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n))));

b = b(r);
n = n(c);

systems(1).burgers=b;
systems(1).plane=n;

burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;


for i=1:size(b,1)    
    t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS);
    t.dispStyle='uvw';
    t=round(t);
    line(i+Ntypes)=t.normalize;
end



Ntypes=size(burgers,2);
edge=1:Ntypes;


%**************************
%Screw Dislocations 
%**************************
b=Miller(1,1,0,CS,'uvw');
[b,c] = symmetrise(b);

systems(2).burgers=b;
systems(2).plane='screw';


burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion;

line(Ntypes+1:Ntypes+size(b,1))=b.normalize;

Ntypes=size(burgers,2);
screw=(edge(end)+1):Ntypes;



%Set up A matrix and other problem parameters

%assign to temporary variables
bs=burgers;
ls=line;
%Calculate A matrix
A = -1*[bs.x.*ls.x-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.x.*bs.y; ls.y.*bs.x; ls.y.*bs.y-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.z.*bs.x; ls.z.*bs.y];

f(edge)=1;
f(screw)=1-poisson;  

lb=zeros(Ntypes,1);

systems(1).indices=edge;
systems(2).indices=screw;

end