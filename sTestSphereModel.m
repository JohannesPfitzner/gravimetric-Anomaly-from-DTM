clear all; %#ok
close all;
clc;

%% Intialization

handle = cHandleData;

handle.density = 2700;
handle.resFactor = 1;
handle.pathData = '../Data/';
handle.pathImage = '../Images/';

% Observation points

handle.xObs = -5000:10:5000;
handle.yObs = zeros(size(handle.xObs));
handle.zObs = zeros(size(handle.xObs));

matObs = [handle.xObs' handle.yObs' handle.zObs'];

% Sphere parameters

r = 500.0;
center = [0,0,1000];
nPhi = 40;
nThetaHalf = 20;

% Create sphere

PHI        = linspace(0, 2*pi, nPhi);
THETA_top  = linspace(0, 0.5*pi, nThetaHalf);
THETA_bot  = linspace(0.5*pi, pi, nThetaHalf);

[phi_top, theta_top] = meshgrid(PHI, THETA_top);
[phi_bot, theta_bot] = meshgrid(PHI, THETA_bot);

x_top = r*sin(theta_top).*cos(phi_top) + center(1);
y_top = r*sin(theta_top).*sin(phi_top) + center(2);
z_top = r*cos(theta_top) + center(3);

x_bot = r*sin(theta_bot).*cos(phi_bot) + center(1);
y_bot = r*sin(theta_bot).*sin(phi_bot) + center(2);
z_bot = r*cos(theta_bot) + center(3);

handle.xSurf = [x_top(:); x_bot(:)];
handle.ySurf = [y_top(:); y_bot(:)];
handle.zSurf = [z_top(:); z_bot(:)];

%% Calculate analytical solution

gx = zeros(length(handle.xObs),1);
gy = gx;
gz = gx;

for i = 1:length(handle.xObs)
    [gx(i), gy(i), gz(i)] = fGravSphere(r,handle.density,center,[handle.xObs(i) handle.yObs(i) handle.zObs(i)]);
end

gx = gx*1e2;
gy = gy*1e2;
gz = gz*1e2;

%% Using Magranaso

% Triangulation

trianglesTop = delaunay(x_top,y_top);
trianglesBot = delaunay(x_bot,y_bot)+(nPhi*nThetaHalf);
trianglesBot_old = trianglesBot;
trianglesBot = [trianglesBot(:,2) trianglesBot(:,1) trianglesBot(:,3)];

handle.triangles = [trianglesTop; trianglesBot];
handle.points = [handle.xSurf, handle.ySurf, handle.zSurf];

%%
% Compute grav anomaly

fGravAnomaly(handle);

handle.Gx = handle.Gx * 1e3;
handle.Gy = handle.Gy * 1e3;
handle.Gz = handle.Gz * 1e3;

%% Using GBOX

% Create cube which has the same volume and center as the sphere

k = nthroot(4/3*pi,3) * r;  % Edgelength

corner_1 = center + [k/2 k/2 k/2];
corner_2 = center + [-k/2 -k/2 -k/2];

% Compute grav anomaly

gzGBOX = zeros(length(handle.xObs),1);

for i = 1:length(handle.xObs)
    gzGBOX(i) = gbox(handle.xObs(i),handle.yObs(i),handle.zObs(i),...
                           corner_2(1),corner_2(2),corner_2(3),...
                           corner_1(1),corner_1(2),corner_1(3),...
                           handle.density);
end




%% Using IGMAS (just write TSRUF)

points = [(1:length(handle.points(:,1)))' handle.points];

%WriteTSurf(handle.pathData,'TsurfSphereModel',points,handle.triangles);

dataIGMAS = csvread([handle.pathData 'calculatedSphereNegDensity.csv'],1,0);

gzIGMAS = dataIGMAS(:,4);

%% Calculate volume

VSphere = 3/4 * pi * r^3;

VTriangulation = 0.0;

for i = 1:size(handle.triangles,1)
    v1 = handle.points(handle.triangles(i,1)) - center;
    v2 = handle.points(handle.triangles(i,2)) - center;
    v3 = handle.points(handle.triangles(i,3)) - center;
    VTriangulation = VTriangulation + 1/6*dot(v1,cross(v2,v3));
end

%% Plot

% Triangulation
figTriangulation = figure;
trisurf(handle.triangles,handle.xSurf,handle.ySurf,handle.zSurf);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
zt = get(gca, 'ZTick');
set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt));
axis equal;

% % Analytical solution
% figAnalytical = figure;
% plot(handle.xObs,gx);
% hold on;
% plot(handle.xObs,gy);
% plot(handle.xObs,gz);
% hold off;
% legend('g_x','g_y','g_z','Location','NorthEast');
% ylabel('g [mGal]');
% xlabel('x [m]');
% grid on;
% set(gca,'GridLineStyle','--');
% 
% % Magranaso
% figMagranaso = figure;
% plot(handle.xObs,handle.Gx);
% hold on;
% plot(handle.xObs,handle.Gy);
% plot(handle.xObs,handle.Gz);
% hold off;
% legend('g_x','g_y','g_z','Location','NorthEast');
% ylabel('g [mGal]');
% xlabel('x [m]');
% grid on;
% set(gca,'GridLineStyle','--');

% Comparison
figComparison = figure;
plot(handle.xObs,gz);
hold on;
plot(handle.xObs,handle.Gz);
plot(handle.xObs,gzGBOX');
plot(handle.xObs,gzIGMAS');
hold off;
legend('g_z analytisch','g_z Magranaso','g_z GBOX','g_z IGMAS+','Location','NorthEast');
ylabel('g_z [mGal]');
xlabel('x [m]');
grid on;
set(gca,'GridLineStyle','--');

% Difference
errMagranaso = abs((gz - handle.Gz')./(gz))*100;
errGBOX = abs((gz - gzGBOX)./(gz))*100;
errIGMAS = abs((gz - gzIGMAS)./(gz))*100;

figDifference = figure;
plot(handle.xObs,errMagranaso);
hold on;
plot(handle.xObs,errGBOX);
plot(handle.xObs,errIGMAS);
hold off;
legend('Error_{Magranaso}','Error_{GBOX}','Error_{IGMAS+}','Location','NorthEast');
ylabel('Rel. Error [%]');
xlabel('x [m]');
grid on;
set(gca,'GridLineStyle','--');

%%
% Save figure to '.png'

%saveas(figTriangulation, [handle.pathImage 'TriangulationSphere.png']);
%saveas(figComparison, [handle.pathImage 'GZSphereComparison.png']);
%saveas(figDifference, [handle.pathImage 'GZSphereError.png']);

%%
% Save observation points to '.csv'

tmp = [handle.xObs' handle.yObs' handle.zObs'];

csvwrite([handle.pathData 'obsSphere.csv'],tmp);

clear tmp;