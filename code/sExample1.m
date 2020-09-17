%%  sExample1.m
%
%   @author     Johannes Pfitzner

% Gravitational constant
gamma = 6.67430*1e-11;

% Normal gravity gradient
gGamma = 0.3086;

% Resolution factor to coarsen the surface grid
resFactor = 10.0;

% Average density [kg/m^3] of the mountain body
density = 2700;

% Observation points
xObs = [384611.4074; 384611.4529; 384611.0635; 384611.4567];
yObs = [5643139.0488; 5643139.4324; 5643139.3648; 5643139.1835]; 
zObs = [428.750000; 362.230000; 327.820000; 281.240000];

% Measured 
gObs = [5617.996; 5625.538; 5629.011; 5633.303];

% Height difference between the observation points
deltaH = abs([zObs(4)-zObs(3); ...
           zObs(3)-zObs(2); ...
           zObs(2)-zObs(1); ...
           zObs(4)-zObs(1)]);

% Load surface grid (from an Esri ASCII raster format-file '.asc')
[xSurf, ySurf, zSurf] = fReadGridFromASC(...
    '..\data\Freiberg_dgm_10m_spac.asc', resFactor);

%%  Calculate topographic reduction and layer density using 'Magranaso' (triangles)

tic
% Triangulate the mountain body
[triangles, points] = fTriangulateFromSurface(xSurf(:),ySurf(:),zSurf(:));

% Calculate topographic reduction
gzMagranaso = fTopographicReductionMagranaso(xObs,yObs,zObs, ...
                                            triangles, points,density);
timeMagranaso = toc;
                                         
% Calculate layer density
gzRedMagranaso = gObs + gzMagranaso;

deltaGMagranaso = [gzRedMagranaso(4)-gzRedMagranaso(3); ...
                gzRedMagranaso(3)-gzRedMagranaso(2); ...
                gzRedMagranaso(2)-gzRedMagranaso(1); ...
                gzRedMagranaso(4)-gzRedMagranaso(1)];    

for i = 1:length(deltaGMagranaso)
    densLayerMagranaso(i) = 1/(4*pi*1e5*gamma)*(gGamma - ...
        deltaGMagranaso(i)./deltaH(i)); %#ok
end

%%  Calculate topographic reduction and layer density using 'GBOX'(pillars)

% Calculate topographic reduction
tic
gzGBOX = fTopographicReductionGBOX(xObs,yObs,zObs, ...
                                  xSurf,ySurf,zSurf,density);
timeGBOX = toc;
                              
% Calculate layer density
gzRedGBOX = gObs + gzGBOX;

deltaGGBOX = [gzRedGBOX(4)-gzRedGBOX(3); ...
                gzRedGBOX(3)-gzRedGBOX(2); ...
                gzRedGBOX(2)-gzRedGBOX(1); ...
                gzRedGBOX(4)-gzRedGBOX(1)];    

for i = 1:length(deltaGGBOX)
    densLayerGBOX(i) = 1/(4*pi*1e5*gamma)*(gGamma - deltaGGBOX(i)./ ...
        deltaH(i)); %#ok
end
