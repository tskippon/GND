function curves = CalcCurvature(ebsd)

    try
        grainId = ebsd.prop.grainId;
        x = ebsd.x;
        y = ebsd.y;
        euls = ebsd.rotations.Euler;
        eul1 = euls(1:end,1);
        eul2 = euls(1:end,2);
        eul3 = euls(1:end,3);
        phase = ebsd.phase;
        dx = max(ebsd.unitCell(:,1))-min(ebsd.unitCell(:,1));
        dy = max(ebsd.unitCell(:,2))-min(ebsd.unitCell(:,2));

        
        
        symType=zeros(max(ebsd.phase),1);
        for i=unique(ebsd.phase(ebsd.phase>0))'
            i
            if(strcmp(ebsd(ebsd.phase==i).CS.lattice,'hexagonal'))
                symType(i)=1;
            elseif(strcmp(ebsd(ebsd.phase==i).CS.lattice,'cubic'))
                symType(i)=2;
            end
        end

                
        [curves] = cpp_curvatures(x,y,eul1,eul2,eul3,grainId,phase,dx,dy,symType);
       
    catch
        error('No grainId stored in the EBSD variable. \n%s\n\n%s\n',...
          'Use the following command to store the grainId within the EBSD data',...
          '[grains,ebsd.grainId] = calcGrains(ebsd)')
    end
    
end