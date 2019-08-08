 function [ f, A, lb, systems] = doBurgersBCC(ebsd,phaseNum,poisson) 
 %% Preamble^M 
 %this code outputs array f (containing energies for all dislocation types) 
 %array A (containing Burgers and Line vector information),  
 %array lb (containing lower bounds for density of each dislocation type -i.e. zero),  
 %Structure systems (containing the possibleslip planes and directions and what family they belong to) 
 Ntypes=0; 
 CS=ebsd(ebsd.phase==phaseNum).CS; 
 
 
 
 
 %%Get units of ebsd coordinates and set burgers vector unit conversion^M 
 %%appropriately.(units of Miller type objects are in Angstroms by default,^M 
 %%so need to convert from that).^M 
 
 
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
 b=Miller(1,1,1,CS,'uvw'); 
 n=Miller(1,1,0,CS,'hkl'); 
 
 
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
 first=1:Ntypes; 
 
 
 %************************** 
 %Second Slip System (Edge) 
 %************************** 
 b=Miller(1,1,1,CS,'uvw'); 
 n=Miller(1,1,2,CS,'hkl'); 
 
 
 
 
 [b,c] = symmetrise(b); 
 [n,c] = symmetrise(n,'antipodal'); 
 
 
 [r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n)))); 
 
 
 b = b(r); 
 n = n(c); 
 
 
 systems(2).burgers=b; 
 systems(2).plane=n; 
 
 
 
 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 
 
 
 
 for i=1:size(b,1)     
     t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS); 
     t.dispStyle='uvw'; 
     t=round(t); 
     line(i+Ntypes)=t.normalize; 
 end 
 
 
 
 
 Ntypes=size(burgers,2); 
 secnd=(first(end)+1):Ntypes; 
 
 
 %************************** 
 %Third Slip System (Edge) 
 %************************** 
 b=Miller(1,1,1,CS,'uvw'); 
 n=Miller(1,2,3,CS,'hkl'); 
 
 
 
 
 [b,c] = symmetrise(b); 
 [n,c] = symmetrise(n,'antipodal'); 
 
 
 [r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n)))); 
 
 
 b = b(r); 
 n = n(c); 
 
 
 systems(3).burgers=b; 
 systems(3).plane=n; 
 
 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 
 
 
 
 for i=1:size(b,1)     
     t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS); 
     t.dispStyle='uvw'; 
     t=round(t); 
     line(i+Ntypes)=t.normalize; 
 end 
 
 
 
 
 Ntypes=size(burgers,2); 
 third=(secnd(end)+1):Ntypes; 
 %************************** 
 %Screw Dislocations  
 %************************** 
 b=Miller(1,1,1,CS,'uvw'); 
 [b,c] = symmetrise(b); 
 
 
 systems(4).burgers=b; 
 systems(4).plane='screw'; 
 
 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 line(Ntypes+1:Ntypes+size(b,1))=b.normalize; 
 
 
 
 
 Ntypes=size(burgers,2); 
 screw=(third(end)+1):Ntypes; 
 
 
 
 
 
 
 %assign to temporary variables 
 bs=burgers; 
 ls=line; 
 %Calculate A matrix 
 A = -1*[bs.x.*ls.x-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.x.*bs.y; ls.y.*bs.x; ls.y.*bs.y-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.z.*bs.x; ls.z.*bs.y]; 
 
 
 
 
 
 
 f([first secnd third])=1; 
 f(screw)=1-poisson;  %Poisson ratio for beta Zr =0.45 
 
 
 lb=zeros(Ntypes,1); 
 
 
 
 
 systems(1).indices=first; 
 systems(2).indices=secnd; 
 systems(3).indices=third; 
 systems(4).indices=screw; 
 end 
