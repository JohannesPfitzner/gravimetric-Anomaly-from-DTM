function [xSurf,ySurf,zSurf] = fReadGridFromASC(filename, resFactor)
%fReadGridFromASC
%                       Reads the surface grid stored in an Esri ASCII 
%                       raster format-file '.asc'. Coarsen the grid using a
%                       a resolution factor.
%
%                       @input parameter:
%                           filename   ..   filename of the given Esri 
%                                           ASCII raster format-file
%                                           ('filename.asc').
%                           resFactor  ..   factor by which the resolution
%                                           of the x-y grid is coarsened
%
%                       @output parameter:
%                           xSurf,ySurf,zSurf
%                                      ..   surface data 
%
%

    [zSurf, metadata] = arcgridread(filename, 'planar');
    
    nrows = metadata.RasterSize(1);
    ncols = metadata.RasterSize(2);
    
    yMin = metadata.YWorldLimits(1) + metadata.YIntrinsicLimits(1);
    yMax = metadata.YWorldLimits(2) - metadata.YIntrinsicLimits(2) + nrows;
    xMin = metadata.XWorldLimits(1) + metadata.XIntrinsicLimits(1);
    xMax = metadata.XWorldLimits(2) - metadata.XIntrinsicLimits(2) + ncols;
    
    if strcmp(metadata.ColumnsStartFrom(), 'north')
        yy = yMax:-metadata.CellExtentInWorldY:yMin;
    elseif strcmp(metadata.ColumnsStartFrom(), 'south')
        yy = yMin:metadata.CellExtentInWorldY:yMax;
    end
    
    if strcmp(metadata.RowsStartFrom(), 'east')
        xx = xMax:-metadata.CellExtentInWorldX:xMin;
    elseif strcmp(metadata.RowsStartFrom(), 'west')
        xx = xMin:metadata.CellExtentInWorldX:xMax;
    end
    
    %xSpacing = metadata.CellExtentInWorldX * resFactor;
    %ySpacing = metadata.CellExtentInWorldY * resFactor;
    
    % Coarsen the grid
    xx = xx(1:resFactor:end);
    yy = yy(1:resFactor:end);
    zSurf = zSurf(1:resFactor:end,1:resFactor:end);
    
    %nRows = nrows/resFactor;
    %nCols = ncols/resFactor;
    
    [xSurf, ySurf] = meshgrid(xx,yy);
end

