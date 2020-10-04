%% sBlock.m
%   Calculate the gravitational force of a cuboid shaped mountain block
%   with the same density and volume as the mountain body described by the
%   surface grid. The floor area is the same.

% Gravitational constant
gamma = 6.67430*1e-11;

% Normal gravity gradient
gGamma = 0.3086;

% Density of the mountain body in kg/m^3
density = 2700; 

% Resolution factor to coarsen the surface grid
resFactor = 1.0;

% Load surface grid (from an Esri ASCII raster format-file '.asc')
[xSurf, ySurf, zSurf] = fReadGridFromASC('..\data\Freiberg_dgm_10m_spac.asc', resFactor);

% Observation points
xObs = [384611.4074; 384611.4529; 384611.0635; 384611.4567];
yObs = [5643139.0488; 5643139.4324; 5643139.3648; 5643139.1835]; 
zObs = [428.750000; 362.230000; 327.820000; 281.240000];

% Measured 
gObs = [5617.996; 5625.538; 5629.011; 5633.303];

deltaH = abs([zObs(4)-zObs(3); ... %P0 - P1
              zObs(3)-zObs(2); ... %P1 - P2
              zObs(2)-zObs(1); ... %P2 - P3
              zObs(4)-zObs(1)]);   %P0 - P3   

%% Calculate der volume of the mountain body

xMin = min(xSurf,[],'all');
xMax = max(xSurf,[],'all');
yMin = min(ySurf,[],'all');
yMax = max(ySurf,[],'all');
zMin = min(zSurf,[],'all');

z = zSurf - zMin;

% Volume of the mountain body defined by the surface grid with 10mx10m 
% cells scaled by resolution factor 
V = (10*resFactor)*(10*resFactor)*sum(z,'all');

% Width and length of the mountain body
xlength = xMax - xMin;
ylength = xMax - xMin;

% Basal area of the mountain body
A = xlength * ylength;

% Height of a cuboid with the same basal area and volume
hCuboid = V/A;

%% Compute the vertical gravitational force of the cuboid

offset = 0.0001;

top = hCuboid + zMin;

gz = zeros(length(zObs),1);

for i = 1:length(zObs)

    if(zObs(i)<zMin)    % Observer point below cuboid
        gz(i) = -gbox(xObs(i),yObs(i),0.0, ...
            xMin,yMin,zMin-zObs(i), ...
            xMax,yMax,(top-zObs(i)), ...
            density);
    elseif(zObs(i)>top) % Observer point above cuboid
        gz(i) = gbox(xObs(i),yObs(i),0.0, ...
            xMin,yMin,(zObs(i)-top), ...
            xMax,yMax,(zObs(i)-zMin), ...
            density);
    else                % Observer point inside cuboid
        % Cuboid below observer point
        gz(i) = +gbox(xObs(i),yObs(i),0.0, ...
            xMin,yMin,offset, ...
            xMax,yMax,(zObs(i)-zMin), ...
            density);
        % Cuboid above observer point
        gz(i) = gz(i) - gbox(xObs(i),yObs(i),0.0, ...
            xMin,yMin,offset, ...
            xMax,yMax,(top-zObs(i)), ...
            density);
    end    
end

%% Calculate layer density

gzRed = gObs + gz;

deltaG = [gzRed(4)-gzRed(3); ...
          gzRed(3)-gzRed(2); ...
          gzRed(2)-gzRed(1); ...
          gzRed(4)-gzRed(1)];    

for i = 1:length(deltaG)
    densLayerGBOX(i) = 1/(4*pi*1e5*gamma)*(gGamma - deltaG(i)./deltaH(i)); %#ok
end
