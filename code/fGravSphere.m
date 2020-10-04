function [gx, gy, gz] = fGravSphere(radius, density, center,...
                                    matObserver)
%FGRAVSPHERE Summary of this function goes here
%   Detailed explanation goes here

    gamma = 6.6742*1e-11;
    m = 4/3 * pi * density * (radius^3);
    r = matObserver - center;
    if(norm(r) > radius)
        g = gamma.*m.*r./(norm(r).^3);
    else
        g = 4/3 * pi * gamma * density * r;
    end
    g(isnan(g)) = 0;
    
    % Convert [m/s^2] to [mGal]
    gx = 1e5*g(:,1);
    gy = 1e5*g(:,2);
    gz = 1e5*g(:,3);
end

