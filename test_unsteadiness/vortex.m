function [u, v, w] = vortex(x, y, z, x1, y1, z1, x2, y2, z2, gamma)

    R_cut = 1e-10;

    %% calculation of R1 x R2
    
    R1R2x = (y - y1) * (z - z2) - (z - z1) * (y - y2);
    R1R2y = -((x - x1) * (z - z2) - (z - z1) * (x - x2));
    R1R2z = (x - x1) * (y - y2) - (y - y1) * (x - x2);
    
    % calculation of (R1 x R2)^2
    square = R1R2x * R1R2x + R1R2y * R1R2y + R1R2z * R1R2z;
    
    %% calculation of R0(R1/R(R1) - R2/R(R2))
    R1 = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));
    R2 = sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2) + (z - z2) * (z - z2));

    if R1 < R_cut || R2 < R_cut || square < R_cut
        u = 0;
        v = 0;
        w = 0;
    else
        R0R1 = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) + (z2 - z1) * (z - z1);
        R0R2 = (x2 - x1) * (x - x2) + (y2 - y1) * (y - y2) + (z2 - z1) * (z - z2);
        coef = gamma / (4 * pi * square) * (R0R1 / R1 - R0R2 / R2);
        u = R1R2x * coef;
        v = R1R2y * coef;
        w = R1R2z * coef;
    end

end
