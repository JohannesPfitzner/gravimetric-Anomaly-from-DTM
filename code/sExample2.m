%% sExample2.m


% Resolution factor to coarsen the surface grid
resFactor = 10.0;

% Read measurement data
dataTable = xlsread('..\Data\Ablesungen_2019.xlsx');

% Time
[~,~,days,hours,minutes,seconds] = datevec([dataTable(3:7,5), ...
                                         dataTable(9:13,5)]);                                  
                                     
time = days*86400 + hours*3600 + minutes*60 + seconds;
 
% dg [mGal]
dg = [dataTable(3:7,3), ...
     dataTable(9:13,3)];

% Density [kg/m^3]
rho = 2700;
 
% Height of the gravimeter [m]
h = 0.2; 

% Height of the measurement points above NHN
hNHN = [281.24; 327.82; 362.23; 428.75];

% Measurement point coordinates
xObs = [384611.0; 384611.0; 384611.0; 384611.0];
yObs = [5643139.0; 5643139.0; 5643139.0; 5643139.0];

% Depth hanging bench [m]
z1 = 0.0;

% Depth 'Rothschönberger Stollen' [m]
z2 = 229.0;

% Depth [m]
z = [dataTable(10,10); ...
     dataTable(17,10); ...
     dataTable(24,10); ...
     dataTable(31,10)];

% Drift height [m]
H = [dataTable(11,10); ...
     dataTable(18,10); ...
     dataTable(25,10)];
 
% Drift width [m]
B = [dataTable(12,10); ...
     dataTable(19,10); ...
     dataTable(26,10)];

% Abstand zum Stoß [m]
x = [dataTable(13,10); ...
     dataTable(20,10); ...
     dataTable(27,10)];

% Schachtquerschnitt [m]
A = [dataTable(14,10); ...
     dataTable(20,10); ...
     dataTable(28,10); ...
     dataTable(32,10)];

% Abstand Gravimeter zu Schachtmittelpunkt [m]
r = [dataTable(15,10); ...
     dataTable(21,10); ...
     dataTable(29,10); ...
     dataTable(33,10)];

t1 = z - z1;
t2 = z - z2;

%%
% Gangkorrektur bestimmen

tmpGGang(1,:) = dg(1,:) - dg(1,1);
tmpGGang(2,:) = dg(5,:) - dg(1,1);

timeGang(1,:) = time(1,:) - time(1,1);
timeGang(2,:) = time(5,:) - time(1,1);

mGang = mean(tmpGGang(2:end)./timeGang(2:end));

gGangBP = mGang * timeGang;

gGang = dg - mGang * time;

%%
% Mittelwert der Mehrfachmessung

gMean(1) = mean([gGang(1,:), gGang(5,:)],2);
gMean(2:4) = mean(gGang(2:4,:),2);
gMean = gMean';

%%
% Anomalie aus Reduktionen berechnen

% Streckenreduktion
deltaGStr = rectprism(x, h, B, H, rho);
deltaGStr(4) = 0.0;

% Schachtreduktion
deltaGSch = vertical_cylinder(r, t1, t2, rho, A);
deltaGSch(4) = 0.0;

% Hohlraumreduktion
deltaGH = zeros(size(deltaGStr));
deltaGH(3,:) = gbox(0.0, 0.0, 0.0, -2.8, -1.3, 0.3, 10.8, 1.3, 9.0, -rho);
%%
% Geländereduktion
[xSurf, ySurf, zSurf] = fReadGridFromASC( '..\data\Freiberg_dgm_10m_spac.asc', resFactor);

[triangles, points] = fTriangulateFromSurface(xSurf(:), ySurf(:), zSurf(:));

deltaGTop = fTopographicReductionMagranaso(xObs,yObs,zObs, ...
                                             triangles,points,density);

% Anomalie
DeltaG = gMean + deltaGSch + deltaGStr + deltaGH + deltaGTop;

%%
% Blockdichte berechnen 

diffDeltaG = [DeltaG(1)-DeltaG(2); ... %P0 - P1
              DeltaG(2)-DeltaG(3); ... %P1 - P2       
              DeltaG(3)-DeltaG(4); ... %P2 - P3
              DeltaG(1)-DeltaG(4)];    %P0 - P3
    
deltaH = abs([hNHN(1)-hNHN(2); ... %P0 - P1
              hNHN(2)-hNHN(3); ... %P1 - P2
              hNHN(3)-hNHN(4); ... %P2 - P3
              hNHN(1)-hNHN(4)]);   %P0 - P3     
    
gamma = 6.67430*1e-11;  % Gravitationskonstante 
g_gamma = 0.3086;   % Normalschweregradient

for i = 1:length(diffDeltaG)
    rhoBMagra(i) = 1/(4*pi*1e5*gamma)*(g_gamma - diffDeltaG(i)./deltaH(i)); %#ok
end

%% Plot

% figGang = figure;%('units','normalized','outerposition',[0 0 1 1]);
% s = scatter(timeGang(:),tmpGGang(:));
% s.MarkerFaceColor = [0 0.5 0.5];
% hold on;
% plot(timeGang(:),gGangBP(:))
% hold off;
% xlabel('Zeit [h]');
% ylabel('g_{Gang} [mGal]');
% grid on;
% set(gca,'GridLineStyle','--');

% figStr = figure;%('units','normalized','outerposition',[0 0 1 1]);
% s = scatter(hNHN,deltaGStr,'o');
% s.MarkerFaceColor = [0 0.5 0.5];
% xlabel('h_{NHN} [m]');
% ylabel('\Delta g_{Str} [mGal]');
% grid on;
% xticks(hNHN');
% xticklabels({'P_0','P_1','P_2','P_3'});
% set(gca,'GridLineStyle','--');

% figSch = figure;%('units','normalized','outerposition',[0 0 1 1]);
% s = scatter(hNHN,deltaGSch,'o');
% s.MarkerFaceColor = [0 0.5 0.5];
% xlabel('h_{NHN} [m]');
% ylabel('\Delta g_{Sch} [mGal]');
% grid on;
% xticks(hNHN');
% xticklabels({'P_0','P_1','P_2','P_3'});
% set(gca,'GridLineStyle','--');

% figDens = figure;%('units','normalized','outerposition',[0 0 1 1]);
% s = scatter(hNHN,rhoBMagra,'o');
% s.MarkerFaceColor = [0 0.5 0.5];
% xlabel('Gesteinsblock');
% ylabel('\rho_B [kg/m^3]');
% grid on;
% xticks(hNHN');
% xticklabels({'P_0 - P_1','P_1 - P_2','P_2 - P_3','P_0 - P_3'});
% set(gca,'GridLineStyle','--');