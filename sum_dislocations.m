 function [GND] = sum_dislocations(disArray, systems,ebsd) 
 %sum_dislocations collects individual dislocation types and group them into 
 %slip systems for easier plotting/analysis 
 %   Detailed explanation goes here 
 
 
 %disArray = output from GND_auto.m containing a matrix of the dislocation 
 %densities for all dislcation types at each point in an ebsd map 
 
 
 %systems = output from GND_auto.m containing information about the burgers 
 %vectors, slip planes, phases, etc. of the dislocation types in disArray 
 
 
 %GND = a structure containing the information from disArray, arranged 
 %according to the slip systems in systems.  Properties of the structure 
 %include burgers, plane, phase, name, data.  Data contains all the GND 
 %densities for that system, name is the name of the slip system (for HCP 
 %systems), plane and burgers are the slip plane and Burgers vectors, and 
 %phase is the phase # 
 
 
 %For example, to get the slip plane of the first slip system, use 
 %GND(1).plane, and to get the dislocation density of dislocations on that 
 %slip system use GND(1).data 
 
 
 
 
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
