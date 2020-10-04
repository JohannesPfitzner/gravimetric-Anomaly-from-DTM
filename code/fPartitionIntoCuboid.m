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
    offset = 0.0001;
    points = [];
    cuboids = [];

    % Consider oberserver point as zero 
    Z = zSurf-zObs;
    
    zMin = min(Z(:));
    
    % Surface grid spacing
    xSpacing = abs(xSurf(1,2)-xSurf(1,1));
    ySpacing = abs(ySurf(2,1)-ySurf(1,1));
    
    % Check for surface points exactly on observer point level and move 
    % them upwards by a fixed offset
    atObserverLevel = (Z==0);
    Z(atObserverLevel) = Z(atObserverLevel) + offset;
    
    % Find the matrix positions of the surface points above and below
    % observer point level
    aboveObserverPoint = Z>0;
    belowObserverPoint = Z<0;   
    
    % Open figure to plot the different pillars by drawing their two
    % defining points
    figure
    
    % Pillars above observer point level
    if(any(aboveObserverPoint(:)))
        tmpx = xSurf(aboveObserverPoint);
        tmpy = ySurf(aboveObserverPoint);
        tmpz = Z(aboveObserverPoint);
        
        pointsTopUp = [tmpx(:), tmpy(:), tmpz(:)];
        pointsTopDown = [tmpx(:) + xSpacing, ...
                         tmpy(:) + ySpacing, ...
                         zeros(length(tmpz(:)),1) + offset];  

        % If observer point is below the deepest surface point
        if(zMin>0.0)
            pointsTopDown(:,3) = pointsTopDown(:,3) + zMin;
        end
                     
        points = [pointsTopUp; pointsTopDown];
        cuboids = [(1:length(tmpx(:)))', ...
                   (1:length(tmpx(:)))' + size(pointsTopUp,1)];
        
        % Plot
        plot3([pointsTopUp(:,1);pointsTopDown(:,1)],...
        [pointsTopUp(:,2);pointsTopDown(:,2)],...
        [pointsTopUp(:,3);pointsTopDown(:,3)],'b.');
        hold on;
               
        % Pillars below the top pillars (when observer point is between the
        % lowest and highest surface level)
        if(any(belowObserverPoint(:)))
            tmpz1 = zeros(size(tmpx));
            tmpz2 = zeros(size(tmpx));
            tmpz1(:) = -offset;
            tmpz2(:) = zMin;

            points1 = [tmpx(:), tmpy(:), tmpz1(:)];
            points2 = [tmpx(:) + xSpacing, ...
                       tmpy(:) + ySpacing, ...
                       tmpz2(:)];  

            cuboids = [cuboids;
                       (1:length(tmpx(:)))' + size(points,1), ...
                       (1:length(tmpx(:)))' + size(points,1) + size(points1,1)];
            points = [points; points1; points2];
            
            % Plot
            plot3([points1(:,1);points2(:,1)],...
            [points1(:,2);points2(:,2)],...
            [points1(:,3);points2(:,3)],'g.');
        end
    end
    
    % Pillars below observer point level
    if(any(belowObserverPoint(:)))
        tmpxx = xSurf(belowObserverPoint);
        tmpyy = ySurf(belowObserverPoint);
        tmpzz = Z(belowObserverPoint);
        
        pointsBotUp = [tmpxx(:), tmpyy(:), tmpzz(:)];
        pointsBotDown = [tmpxx(:) + xSpacing, ...
                         tmpyy(:) + ySpacing, ...
                         zeros(length(tmpzz(:)),1) + zMin];  

        cuboids = [cuboids;
                   (1:length(tmpxx(:)))' + size(points,1), ...
                   (1:length(tmpxx(:)))' + size(points,1) + size(pointsBotUp,1)];
        points = [points; pointsBotUp; pointsBotDown];
        
        % Plot
        plot3([pointsBotUp(:,1);pointsBotDown(:,1)],...
        [pointsBotUp(:,2);pointsBotDown(:,2)],...
        [pointsBotUp(:,3);pointsBotDown(:,3)],'r.');      
    end
    
    % Plot options
    if(all(aboveObserverPoint(:))) 
        hold off;
        legend('Top');
    elseif(all(belowObserverPoint(:)))
        legend('Bot');        
    else
        hold off;
        legend('Top','TopBelow','Bot');
    end
    title(['Pillars for observer point at z = ' int2str(zObs) 'm']);
end

