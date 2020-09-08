%%  Init

% Observation points
xObs = [384611.4074; 384611.4529; 384611.0635; 384611.4567];
yObs = [5643139.0488; 5643139.4324; 5643139.3648; 5643139.1835]; 
zObs = [428.750000; 362.230000; 327.820000; 281.240000]; 

% Filename of the surface data
filename = 'Freiberg_dgm_10m_spac.asc';

% Average density [kg/m^3] of the mountain body
density = 2700;

% Resolution factor to coarsen the surface grid
resFactor = 10.0;

% Load surface grid (from an Esri ASCII raster format-file '.asc')
[xSurf, ySurf, zSurf] = fReadGridFromASC(filename, resFactor);

%%  Calculate topographic reduction using 'Magranaso' (triangles)

% Triangulate mountain body using surface data
[triangles, points] = fTriangulateFromSurface(xSurf, ySurf, zSurf);

% Compute topographic reduction
gzMagranaso = fTopographicReductionMagranaso(xObs,yObs,zObs, ...
                                             triangles,points,density);

%%  Calculate topographic reduction using 'GBOX' (pillars)

gzGBOX = fTopographicReductionGBOX(xObs,yObs,zObs, ...
                                   xSurf,ySurf,zSurf,density);