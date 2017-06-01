function [disArray, systems]= GND_auto( ebsd, nthreads,poisson,cubicType )
%Calculate GND densities for ebsd map using parallel procesing.
%ebsd= ebsd data with grain reconstruction complete (reccomended that the data is smoothed first)
%nthreads= number of threads you wish to use for parallel computing
%poisson= an array containing the poisson's ratios for each phase (if this
%is missing you will be prompted to enter the information)
%disArray= array containing GND data
%systems= structure containing info about how different dislocation types
%are stored in disArray


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


%Calculate lattice curvatures using C language code for speed.
disp('Calculating Curvatures...');
curve = CalcCurvature(ebsd);
disp('Curvature calculations complete!');

    
tempPoisson=poisson;
%automated setup of dislocation types
for i=unique(ebsd.phase(ebsd.phase>0))'
    if(length(unique(ebsd.phase(ebsd.phase>0)))>poissonDefined)
        prompt=sprintf('Number of defined Poissons Ratio less than number of phases in ebsd data.  Enter Poisson Ratio for phase %i (%s) now.\n (To avoid this dialog send an array containing the values for all phases as an input when calling GND code)',i,ebsd(ebsd.phase==i).mineral);
        name = 'Poissons ratio:';
        defaultans = {'0.30'};
        input = inputdlg(prompt,name,[1 40],defaultans);
        poisson(i)=str2double(input{:});
    else
        poisson(i)=tempPoisson(find(unique(ebsd.phase(ebsd.phase>0))'==i));
    end
    if(strcmp(ebsd(ebsd.phase==i).CS.lattice,'hexagonal'))
     
        [f{i}, A{i}, lb{i},systems{i}]=doBurgersHCP(ebsd,i,poisson(i));
    end
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



%**************************************************************************
% Minimization
%**************************************************************************
disp('Minimizing dislocation energy...');


singleCoreTime=40/2000; %single core takes 40 seconds for 2000 pointss
speedUp=0.5305*nthreads+0.7107; %Amount of speedup gained from using nthreads cores
if nthreads==1
    speedUp=1;
end
estimate=size(curve,1)*singleCoreTime/speedUp;
fprintf('Estimated time to completion: %s if running on %d cores.\n',datestr(estimate/60/60/24,'HH:MM:SS'),nthreads);
%assignin('base', 'curve', curve);
disArray = zeros(size(curve,1),max(cellfun('length',f)));

alpha1=-curve(:,2);
alpha2=-curve(:,3);

alpha3=-curve(:,4);
alpha4=-curve(:,6);
alpha5=curve(:,1)+curve(:,5);

options=optimoptions('linprog','Algorithm','dual-simplex','Display','off');

tic
phase=ebsd.phase;
fmaxSize=max(cellfun('length',f));

%loop through all points
parfor (j=1:size(curve,1),nthreads)
    %for j=1:size(curve,1) %this line for running in serial
    if(phase(j)~=0)
        x =linprog(f{phase(j)},A{phase(j)},[alpha1(j) alpha2(j) alpha3(j) alpha4(j) alpha5(j)],[],[],lb{phase(j)},[],[],options);
        disArray(j,:) = [x; zeros(fmaxSize-length(f{phase(j)}),1)];
    end
end


disp('Minimization complete!');
actualTime=toc;
estimateError=(estimate-actualTime)/estimate*100;
fprintf('Actual time taken was %s.  Estimate was off by %3.1f%%\n',datestr(actualTime/60/60/24,'HH:MM:SS'),estimateError)
end

