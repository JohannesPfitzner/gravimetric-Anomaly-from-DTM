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
        % observation point z-level
        [cuboids, points] = fPartitionIntoCuboid(xSurf,ySurf,zSurf,zObs(i));
        
        % Calculate the topographic reduction
        sumTop = 0.0;
        sumBot = 0.0;
        for j = 1:size(cuboids,1)
            if(points(cuboids(j,1),3)>0)
                sumTop = sumTop + gbox(xObs(i),yObs(i),0.0, ...
                     points(cuboids(j,1),1), ...
                     points(cuboids(j,1),2), ...
                     points(cuboids(j,1),3), ...
                     points(cuboids(j,2),1), ...
                     points(cuboids(j,2),2), ...
                     points(cuboids(j,2),3), ...
                     density);
            else
                sumBot = sumBot + gbox(xObs(i),yObs(i),0.0, ...
                     points(cuboids(j,1),1), ...
                     points(cuboids(j,1),2), ...
                     -points(cuboids(j,1),3), ...
                     points(cuboids(j,2),1), ...
                     points(cuboids(j,2),2), ...
                     -points(cuboids(j,2),3), ...
                     density);
            end
        end
        Gz(i) = sumTop + sumBot;
        waitbar(i/length(zObs),h);
    end
    close(h);
    
end

% function g = gbox(x0, y0, z0, x1, y1, z1, x2, y2, z2, rho)
% %gbox
% % c  Subroutine GBOX computes the vertical attraction of a 
% % c  rectangular prism.  Sides of prism are parallel to x,y,z axes,
% % c  and z axis is vertical down.  
% % c
% % c  Input parameters:
% % c    Observation point is (x0,y0,z0).  The prism extends from x1
% % c    to x2, from y1 to y2, and from z1 to z2 in the x, y, and z 
% % c    directions, respectively.  Density of prism is rho.  All 
% % c    distance parameters in units of km; rho in units of 
% % c    kg/(m**3). 
% % c
% % c  Output parameters:
% % c    Vertical attraction of gravity, g, in mGal.
% % 
% % Implemented after R. J. Blakeley, Potential theory in gravity and magnetic applications
% %
% 
% f = 6.673e-11;
% isign = [-1.0, 1.0];
% x = zeros(2, 1);
% y = zeros(2, 1);
% z = zeros(2, 1);
% x = [x0, x0] - [x1, x2];
% y = [y0, y0] - [y1, y2];
% z = [z0, z0] - [z1, z2];
% sum = 0.0;
% for i = 1:2
%     for j = 1:2
%         for k = 1:2
%             rijk = sqrt(x(i)^2+y(j)^2+z(k)^2);
%             ijk = isign(i) * isign(j) * isign(k);
%             arg1 = atan2((x(i)*y(j)),(z(k)*rijk));
%             if arg1 < 0.0
%                 arg1 = arg1 + 2.0 * pi;
%             end
%             arg2 = rijk + y(j);
%             arg3 = rijk + x(i);
%             arg2 = log(arg2);
%             arg3 = log(arg3);
%             sum = sum + ijk * (z(k) * arg1 - x(i) * arg2 - y(j) *  arg3);
%         end
%     end
% end
% g = rho * f * sum * 1.0e5;
% end
