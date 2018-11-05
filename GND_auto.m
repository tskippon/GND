function [disArray, systems]= GND_auto( ebsd, nthreads,poisson,cubicType )
%Calculate GND densities for ebsd map using parallel procesing.
%ebsd= ebsd data with grain reconstruction complete (reccomended that the data is smoothed first)
%nthreads= number of threads you wish to use for parallel computing
%poisson= an array containing the poisson's ratios for each phase (if this
%is missing you will be prompted to enter the information)
%disArray= array containing GND data
%systems= structure containing info about how different dislocation types
%are stored in disArray
%cubicType = cell array setting either BCC or FCC for any cubic phases, in
%order.  Non-cubic phases should NOT be specified.  For example, if phase 1 
%is BCC, phase 2 is HCP, and phase 3 is FCC, then cubicType={'BCC' 'FCC'}


%Start keeping track of elapsed time
tic;

%Set up parallel pool
if nthreads>1
    pool=gcp;
end

if nargin==2
    poissonDefined=0;
end
if nargin>2
    poissonDefined=length(poisson);
end
if nargin>3
    numCubics=0;
    cubicTypeDefined=length(cubicType);
    for i=unique(ebsd.phase(ebsd.phase>0))'
        if strcmp(ebsd(ebsd.phase==i).CS.lattice,'cubic')
            numCubics=numCubics+1;
        end
    end
end

%cubic counter used for checking that cubictype is defined for each cubic
%phase.  Should be set at 1 initially. If there are no cubic phases this
%will not be used.
cubicCounter=1;

%Allocate memory for temporary array to hold curvature data
tempcurve=zeros(length(ebsd),6);
tempcurve(:,:)=NaN;

    
tempPoisson=poisson;
%automated setup of dislocation types

%loop through all phases except phase 0, which is non-indexed data
for i=unique(ebsd.phase(ebsd.phase>0))'
    
    %For each phase, gridify the ebsd data and get the x and y gradients
    temp=ebsd(ebsd.phase==i);
     [temp,newId]=gridify(temp);
    
    gx=temp.gradientX;
    gx=gx(temp.phase==i);
    
    gy=temp.gradientY;
    gy=gy(temp.phase==i);
    
    %Put the curvature components into a temporary variable
    tempcurve(temp.phase==i,:)=[gx.x gx.y gx.z gy.x gy.y gy.z];
    
    
    
    %If user didn't input enough Poissson's ratio values, then ask for them
    if(length(unique(ebsd.phase(ebsd.phase>0)))>poissonDefined)
        prompt=sprintf('Number of defined Poissons Ratio less than number of phases in ebsd data.  Enter Poisson Ratio for phase %i (%s) now.\n (To avoid this dialog send an array containing the values for all phases as an input when calling GND code)',i,ebsd(ebsd.phase==i).mineral);
        name = 'Poissons ratio:';
        defaultans = {'0.30'};
        input = inputdlg(prompt,name,[1 40],defaultans);
        poisson(i)=str2double(input{:});
    else
        poisson(i)=tempPoisson(find(unique(ebsd.phase(ebsd.phase>0))'==i));
    end
    
    %If current phase is hexagonal then run HCP Burgers vector setup (see
    %doBurgersHCP.m)
    if(strcmp(ebsd(ebsd.phase==i).CS.lattice,'hexagonal'))
        [f{i}, A{i}, lb{i},systems{i}]=doBurgersHCP(ebsd,i,poisson(i));
    end
    
    %If current phase is cubic, then check CubicTypes and run either BCC or
    %FCC Burgers vector setup (see doBurgersBCC.m and doBurgersFCC.m)
    if(strcmp(ebsd(ebsd.phase==i).CS.pointGroup,'m-3m'))
        if(cubicTypeDefined<cubicCounter)
        question=sprintf('Cubic phase detected for phase %i (%s).  Is this phase FCC or BCC?',i,ebsd(ebsd.phase==i).mineral);
        default='BCC';
        cubicType{cubicCounter} = questdlg(question,'Specify cubic type','FCC','BCC',default);
        end
        
        if(strcmp(cubicType{cubicCounter},'BCC'))
            fprintf('BCC structure selected for phase %i (%s)\n',i,ebsd(ebsd.phase==i).mineral);
            [f{i}, A{i}, lb{i},systems{i}]=doBurgersBCC(ebsd,i,poisson(i));
        elseif(strcmp(cubicType{cubicCounter},'FCC'))
            [f{i}, A{i}, lb{i},systems{i}]=doBurgersFCC(ebsd,i,poisson(i));
        end
        cubicCounter=cubicCounter+1;    
    end
    
end


%This makes sure that the curvature data calculated from the gridified
%ebsd map lines up with the original (non-gridified) map.  This 
%allows the results to be easily plotted on the original ebsd map
[temp,newId]=gridify(ebsd);
indx=temp.id2ind(newId);
curve=tempcurve(indx,:);


%**************************************************************************
% Minimization
%**************************************************************************
disp('Minimizing dislocation energy...');

%initialize variable for holding the dislocation densities of each type of
%dislocation for each point in the map.
disArray = zeros(size(curve,1),max(cellfun('length',f)));

%Define components of curvature tensor
kappa21=-curve(:,2);
kappa31=-curve(:,3);

kappa12=-curve(:,4);
kappa32=-curve(:,6);

kappa11=-curve(:,1);
kappa22=-curve(:,5);

%Set up options for linear programming solver
options=optimoptions('linprog','Algorithm','dual-simplex','Display','off');
phase=ebsd.phase;
fmaxSize=max(cellfun('length',f));

%Do short test run on the first 1000*nthreads points to estimate time
%required for the full dataset.
testrun=1:min(1000*nthreads,length(ebsd));
tic
parfor (j=testrun,nthreads)
    if(phase(j)~=0 && sum(isnan(curve(j,:)))==0)
        x =linprog(f{phase(j)},A{phase(j)},[kappa11(j) kappa12(j) kappa21(j) kappa22(j) kappa31(j) kappa32(j)],[],[],lb{phase(j)},[],[],options);
        disArray(j,:) = [x; zeros(fmaxSize-length(f{phase(j)}),1)];
    end
end
estimate=toc*(length(curve)-length(testrun))/length(testrun);

%Ptrint out time estimate in HH:MM:SS format.
%fprintf('Estimated time to completion: %s.\n',datestr(estimate/60/60/24,'HH:MM:SS'));
dt=datetime('now')+seconds(estimate);
fprintf('Analysis should be complete by: %s system clock time.\n',datestr(dt));



%loop through all points
parfor (j=1:size(curve,1),nthreads)
%for j=1:size(curve,1) %Replace above line with this one for running
%wihtout parallel computing

    %Only perform calculations on indexed phases, don't redo calculations
    %that were done in the test run.
    if(phase(j)~=0 && max(testrun==j)==0 && sum(isnan(curve(j,:)))==0)        
        x =linprog(f{phase(j)},A{phase(j)},[kappa11(j) kappa12(j) kappa21(j) kappa22(j) kappa31(j) kappa32(j)],[],[],lb{phase(j)},[],[],options);
        %Place solved dislocation densities for various dislocation types
        %in disArray
        disArray(j,:) = [x; zeros(fmaxSize-length(f{phase(j)}),1)];
    end
end

%Print out system time that the analysis finished at, and total elapsed
%time.
disp('Minimization complete!');
fprintf('Analysis completed at %s system clock time.\n',datestr(datetime('now')))
fprintf('Total elapsed time was: %s \n',datestr(toc/60/60/24,'HH:MM:SS'))
end

