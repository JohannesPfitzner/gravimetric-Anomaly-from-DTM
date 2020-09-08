function Gz = fTopographicReductionGBOX(xObs,yObs,zObs,xSurf,ySurf,zSurf,density)
%fTopographicReductionGBOX 
%                           asd
    Gz = zeros(length(xObs),1);
    
    for i = 1:length(xObs)
        sum = 0.0;
        
        isEqual = zSurf==xObs(i);
        zSurf(isEqual) = zSurf(isEqual)*0.99999;
        
        isTop = zSurf>xObs(i);
        isBot = ~isTop;
        
        if(all(isTop)||all(isBot))
            % All pillars below or all above
            corner2 = [xSurf+0.5*xSpacing ...
                       ySurf+0.5*ySpacing ...
                       zSurf];
            corner1 = [xSurf-0.5*xSpacing ...
                       ySurf-0.5*ySpacing ...
                       ones(max(size(corner2)),1)*min(zSurf(:))];
                  
            for j = 1:max(size(corner1))   
                if(all(isTop))
                    sum = sum + gbox(xObs(i),yObs(i),zObs(i),...
                                 corner1(j,1),corner1(j,2),corner1(j,3),...
                                 corner2(j,1),corner2(j,2),corner2(j,3),...
                                 density);
                else
                    sum = sum - gbox(xObs(i),yObs(i),zObs(i),...
                                 corner1(j,1),corner1(j,2),corner1(j,3),...
                                 corner2(j,1),corner2(j,2),corner2(j,3),...
                                 density);
                end
            end 
        else
            % Mixed pillars
            
            % Pillars below
            cornerBot2 = [xSurf+0.5*xSpacing ...
                          ySurf+0.5*ySpacing ...
                          zSurf];
            cornerBot2(isTop) = xObs(i)*0.99999;
            
            cornerBot1 = [xSurf-0.5*xSpacing ...
                          ySurf-0.5*ySpacing ...
                          ones(max(size(corner1)),1) * min(zSurf(:))];
            
            for j = 1:max(size(cornerBot2))   
                sum = sum - gbox(xObs(i),yObs(i),zObs(i),...
                                 cornerBot1(j,1),cornerBot1(j,2),cornerBot1(j,3),...
                                 cornerBot2(j,1),cornerBot2(j,2),cornerBot2(j,3),...
                                 density);
            end 
                   
            % Pillars above
            cornerTop2 = [xSurf(isTop)+0.5*xSpacing ...
                          ySurf(isTop)+0.5*ySpacing ...
                          zSurf(isTop)];
            tmpTop = ones(max(size(corner1)),1) * xObs(i)*1.00001; 
            cornerTop1 = [xSurf(isTop)-0.5*xSpacing ...
                          ySurf(isTop)-0.5*ySpacing ...
                          tmpTop];
            
            for j = 1:max(size(cornerTop2))   
                sum = sum + gbox(xObs(i),yObs(i),zObs(i),...
                                 cornerTop1(j,1),cornerTop1(j,2),cornerTop1(j,3),...
                                 cornerTop2(j,1),cornerTop2(j,2),cornerTop2(j,3),...
                                 density);
            end 
                      
                   
        end
        
        Gz(i) = sum;
    end
    
end

