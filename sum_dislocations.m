function [GND] = sum_dislocations(disArray, systems,ebsd)
%sum_dislocatiosn Collect individual dislocation types and group them into
%slip systems for easier plotting/analysis
%   Detailed explanation goes here
GND_counter=1;

for i=1:length(systems)


j=1;
while j<=length(systems{i})
    if(isfield(systems{i},'name'))
        GND(GND_counter).name=systems{i}(j).name;
    end
    GND(GND_counter).data(ebsd.phase==i)=sum(disArray(ebsd.phase==i,systems{i}(j).indices),2);
    GND(GND_counter).data(ebsd.phase==0 | ebsd.phase >i)=NaN;
    GND(GND_counter).burgers=systems{i}(j).burgers;
    GND(GND_counter).plane=systems{i}(j).plane;   
    GND(GND_counter).phase=i;
    j=j+1;
    GND_counter=GND_counter+1;
end

GND(GND_counter).name='total';
GND(GND_counter).data(ebsd.phase==i)=sum(disArray(ebsd.phase==i,:),2);
GND(GND_counter).data(ebsd.phase~=i)=NaN;
GND(GND_counter).burgers=[];
GND(GND_counter).plane=[];
GND(GND_counter).phase=i;
GND_counter=GND_counter+1;
end

end

