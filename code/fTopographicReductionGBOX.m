function Gz = fTopographicReductionGBOX(xObs,yObs,zObs,xSurf,ySurf,zSurf,density)
%fTopographicReductionGBOX 
%                           Takes a MxN surface grid (M,N > 1) and an 
%                           observer points z-coordinate to partition the 
%                           body (defined by the surface grid on top and a 
%                           plane, with the the z-level of the lowest
%                           surface height) into pillars above and below 
%                           the observer point z-level.
%   @parameters
%               xObs, yObs, zObs    ..  observation points
%               xSurf, ySurf, zSurf ..  surface grid
%               density             ..  approximate density of the body
%                                       confined by the surface grid (top) 
%                                       and a z-plane (bottom). The planes'
%                                        height is defined by the minimal
%                                       z-value of the surface grid.
%
%   @author     Johannes Pfitzner
%
    Gz = zeros(length(zObs),1);
    
    h=waitbar(0,'Calculating with GBOX...');
    for i = 1:length(zObs)
        % Partition the mountain body into cuboids speratated by
        % observation point height
        [cuboids, points] = fPartitionIntoCuboid(xSurf,ySurf,zSurf,...
            zObs(i));
        
        % Calculate the topographic reduction
        sumTop = 0.0;
        sumBot = 0.0;
        for j = 1:size(cuboids,1)
            % Cuboids above observer point height
            if(points(cuboids(j,1),3)>0)
                sumTop = sumTop + abs(gbox(xObs(i),yObs(i),0.0, ...
                     points(cuboids(j,1),1), ...
                     points(cuboids(j,1),2), ...
                     points(cuboids(j,1),3), ...
                     points(cuboids(j,2),1), ...
                     points(cuboids(j,2),2), ...
                     points(cuboids(j,2),3), ...
                     density));
            % Cuboids below observer point height
            else
                sumBot = sumBot + abs(gbox(xObs(i),yObs(i),0.0, ...
                     points(cuboids(j,1),1), ...
                     points(cuboids(j,1),2), ...
                     -points(cuboids(j,1),3), ...
                     points(cuboids(j,2),1), ...
                     points(cuboids(j,2),2), ...
                     -points(cuboids(j,2),3), ...
                     density));
            end
        end
        Gz(i) = - sumTop + sumBot;
        waitbar(i/length(zObs),h);
    end
    close(h);
    
end