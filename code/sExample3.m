%%  sExample3.m     
%  -- Sphere --

% Observation points

xObs = -5000:10:5000;
yObs = zeros(size(xObs));
zObs = zeros(size(xObs));

% Sphere parameters

density = 2700;
r = 500.0;
center = [0,0,1000];
nPhi = 40;
%nThetaHalf = 20;

% Create sphere

% PHI        = linspace(0, 2*pi, nPhi);
% THETA_top  = linspace(0, 0.5*pi, nThetaHalf);
% THETA_bot  = linspace(0.5*pi, pi, nThetaHalf);
% 
% [phiTop, thetaTop] = meshgrid(PHI, THETA_top);
% [phiBot, thetaBot] = meshgrid(PHI, THETA_bot);
% 
% xTop = r*sin(thetaTop).*cos(phiTop) + center(1);
% yTop = r*sin(thetaTop).*sin(phiTop) + center(2);
% zTop = r*cos(thetaTop) + center(3);
% 
% xBot = r*sin(thetaBot).*cos(phiBot) + center(1);
% yBot = r*sin(thetaBot).*sin(phiBot) + center(2);
% zBot = r*cos(thetaBot) + center(3);
% 
% xSurf = [xTop(:); xBot(:)];
% ySurf = [yTop(:); yBot(:)];
% zSurf = [zTop(:); zBot(:)];

[xSurf,ySurf,zSurf] = sphere(nPhi);

xSurf = r * xSurf(:) + center(1);
ySurf = r * ySurf(:) + center(2);
zSurf = r * zSurf(:) + center(3);

ppoints = [xSurf, ySurf, zSurf];
ppoints = unique(ppoints,'rows');

xTop = ppoints(ppoints(:,3)<=center(3),1);
yTop = ppoints(ppoints(:,3)<=center(3),2);
zTop = ppoints(ppoints(:,3)<=center(3),3);
xBot = ppoints(ppoints(:,3)>=center(3),1);
yBot = ppoints(ppoints(:,3)>=center(3),2);
zBot = ppoints(ppoints(:,3)>=center(3),3);

%% Calculate analytical solution

gx = zeros(length(xObs),1);
gy = gx;
gz = gx;

for i = 1:length(xObs)
    [gx(i), gy(i), gz(i)] = fGravSphere(r,density,center,[xObs(i) yObs(i) zObs(i)]);
end

%% Calculate Magranaso solution

trianglesTop = delaunay(xTop,yTop);
trianglesBot = delaunay(xBot,yBot)+size(xTop,1);
trianglesBot = [trianglesBot(:,2) trianglesBot(:,1) trianglesBot(:,3)];

points = [xTop,yTop,zTop;xBot,yBot,zBot];
triangles = [trianglesTop; trianglesBot];

gzMagranaso = fTopographicReductionMagranaso(xObs,yObs,zObs,triangles,points,density);

%% Calculate GBOX solution

% Create cube which has the same volume and center as the sphere

% Edgelength
k = nthroot(4/3*pi,3) * r;  

corner_1 = center + [k/2 k/2 k/2];
corner_2 = center + [-k/2 -k/2 -k/2];

% Compute grav anomaly

gzGBOX = zeros(length(xObs),1);

for i = 1:length(xObs)
    gzGBOX(i) = gbox(xObs(i),yObs(i),zObs(i),...
                           corner_2(1),corner_2(2),corner_2(3),...
                           corner_1(1),corner_1(2),corner_1(3),...
                           density);
end


%% Calculate volume

VSphere = 3/4 * pi * r^3;

VTriangulation = 0.0;

for i = 1:size(triangles,1)
    v1 = points(triangles(i,1),:) - center;
    v2 = points(triangles(i,2),:) - center;
    v3 = points(triangles(i,3),:) - center;
    VTriangulation = VTriangulation + abs(1/6*dot(cross(v2,v3),v1));
end