function fTriangulateSurf(handle)
%FTRIANGULATEDGM Summary of this function goes here
%   Detailed explanation goes here
    ext_x = [min(handle.xSurf) max(handle.xSurf)];
    ext_y = [min(handle.ySurf) max(handle.ySurf)];
    %ext_z = [min(handle.zSurf) max(handle.zSurf)];
    ext_z = [-1000 max(handle.zSurf)];

    
    offsetz = 1;
    
    % Top
    triangles1 = delaunay(handle.xSurf,handle.ySurf);
    points1 = [[1:length(handle.xSurf)]' handle.xSurf handle.ySurf ...
        handle.zSurf];

    % Floor
    bottomPoints = points1;
    bottomPoints(:,4) = ext_z(1)-offsetz;

    triangles2 = delaunay(bottomPoints(:,2),bottomPoints(:,3));

    % First side
    tmp = points1(points1(:,2)==ext_x(1),2:4);
    tmp2 = bottomPoints(bottomPoints(:,2)==ext_x(1),2:4);
    bs_points1 = [tmp;tmp2];

    triangles3 = delaunay(bs_points1(:,2),...
        [ones(length(tmp),1)*mean(tmp(:,3));...
        tmp2(:,3)]);

    % Second side
    tmp = points1(points1(:,2)==ext_x(2),2:4);
    tmp2 = bottomPoints(bottomPoints(:,2)==ext_x(2),2:4);
    bs_points2 = [tmp;tmp2];

    triangles4 = delaunay(bs_points2(:,2),...
        [ones(length(tmp),1)*mean(tmp(:,3));...
        tmp2(:,3)]);

    % Third side
    tmp = points1(points1(:,3)==ext_y(1),2:4);
    tmp2 = bottomPoints(bottomPoints(:,3)==ext_y(1),2:4);
    bs_points3 = [tmp;tmp2];

    triangles5 = delaunay(bs_points3(:,1),...
        [ones(length(tmp),1)*mean(tmp(:,3));...
        tmp2(:,3)]);

    % Fourth side
    tmp = points1(points1(:,3)==ext_y(2),2:4);
    tmp2 = bottomPoints(bottomPoints(:,3)==ext_y(2),2:4);
    bs_points4 = [tmp;tmp2];

    triangles6 = delaunay(bs_points4(:,1),...
        [ones(length(tmp),1)*mean(tmp(:,3));...
        tmp2(:,3)]);

    % Swap normal vectors
    triangles2 = triangles2(:,[1 3 2]);
    triangles3 = triangles3(:,[1 3 2]);
    triangles6 = triangles6(:,[1 3 2]);

    % Merge and fix notation
    points = points1;
    triangles = triangles1; 
    
    bottomPoints(:,1) = bottomPoints(:,1)+length(points);
    triangles2 = triangles2+length(points);
    points = [points;bottomPoints];
    triangles = [triangles;triangles2];

    bs_points1 = [[1:length(bs_points1)]'+length(points) bs_points1];
    triangles3 = triangles3+length(points);
    points = [points;bs_points1];
    triangles = [triangles;triangles3];

    bs_points2 = [[1:length(bs_points2)]'+length(points) bs_points2];
    triangles4 = triangles4+length(points);
    points = [points;bs_points2];
    triangles = [triangles;triangles4];

    bs_points3 = [[1:length(bs_points3)]'+length(points) bs_points3];
    triangles5 = triangles5+length(points);
    points = [points;bs_points3];
    triangles = [triangles;triangles5];

    bs_points4 = [[1:length(bs_points4)]'+length(points) bs_points4];
    triangles6 = triangles6+length(points);
    points = [points;bs_points4];
    triangles = [triangles;triangles6];
    
    handle.triangles = triangles;
    handle.points = points(:,2:end);
    
    WriteTSurf(handle.pathData,['topoModelRZ' int2str(handle.resFactor...
        )],points,triangles);
end
