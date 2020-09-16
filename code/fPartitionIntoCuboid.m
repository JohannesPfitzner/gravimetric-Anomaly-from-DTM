function [cuboids,points] = fPartitionIntoCuboid(xSurf,ySurf,zSurf,zObs)
%fPartitionIntoCuboid 
%                           Takes a MxN surface grid (M,N > 1) and an 
%                           observer points z-coordinate to partition the 
%                           body (defined by the surface grid on top and a 
%                           z-plane, with the the z-level of the lowest
%                           surface height) into pillars above and below 
%                           the observer point z-level.
%   @input parameters
%               xSurf, ySurf, zSurf ..  surface grid
%               zObs                ..  z-coordinate of observer point
%
%   @output parameters
%               cuboids             ..  Lx2 Matrix, where L is the number 
%                                       of cuboids. Each row contains the 
%                                       references to the position in the 
%                                       'points' list the two opposing 
%                                       edges of a cuboid have. (P1,P2)
%               points              ..  Kx3 Matrix, containing all the
%                                       cuboid points, where K is the count
%                                       of points. (x,y,z)
%
%   @author     Johannes Pfitzner
%

    % Consider oberserver point level the zero level
    Z = zSurf-zObs;
    zMin = min(Z,[],'all');
    
    % Surface grid spacing
    xSpacing = abs(xSurf(1,2)-xSurf(1,1));
    ySpacing = abs(ySurf(2,1)-ySurf(1,1));
    
    % Check for surface points exactly on observer point level   
    [row,col] = find(~Z);
    if(~isempty(row))
       offsetFactor = 1.0001;
       Z(row,col) = offsetFactor*Z(row,col);
    end
    
    offset = 0.0001;
    points = [];
    cuboids = [];
    
    % Pillars above observer point level
    if(any(Z>0,'all'))
        tmpx = xSurf(Z>0);
        tmpy = ySurf(Z>0);
        tmpz = Z(Z>0);
        
        pointsTopUp = [tmpx(:), tmpy(:), tmpz(:)];
        pointsTopDown = [tmpx(:) + xSpacing, ...
                         tmpy(:) + ySpacing, ...
                         zeros(length(tmpz(:)),1) + offset];  

        points = [pointsTopUp; pointsTopDown];
        cuboids = [(1:length(tmpx(:)))', ...
                   (1:length(tmpx(:)))' + size(pointsTopUp,1)];
    end
    
    % Pillars below observer point level
    if(any(Z<0,'all'))
        tmpzz = zeros(size(Z));
        tmpzz(Z<0) = Z(Z<0);
        tmpzz(Z>0) = -offset;
        
        pointsBotUp = [xSurf(:), ySurf(:), tmpzz(:)];
        pointsBotDown = [xSurf(:) + xSpacing, ...
                         ySurf(:) + ySpacing, ...
                         zeros(length(tmpzz(:)),1) + zMin];  

        cuboids = [cuboids;
                   (1:length(Z(:)))' + size(points,1), ...
                   (1:length(Z(:)))' + size(points,1) + size(pointsBotUp,1)];
        points = [points; pointsBotUp; pointsBotDown];
    end
    
end

