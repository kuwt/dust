function [A, B] = influence(xc, yc, zc, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
    
    % doublet (a) source (b) influence at point (xc,yc,zc) due to panel (x1,...,z4)

    ep = 1e-6; % panel side cutoff distance

    % distance (R)
    r1 = sqrt((xc - x1)^2 + (yc - y1)^2 + zc^2);
    r2 = sqrt((xc - x2)^2 + (yc - y2)^2 + zc^2);
    r3 = sqrt((xc - x3)^2 + (yc - y3)^2 + zc^2);
    r4 = sqrt((xc - x4)^2 + (yc - y4)^2 + zc^2);

    %panel side
    d1 = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    d2 = sqrt((x3 - x2)^2 + (y3 - y2)^2);
    d3 = sqrt((x4 - x3)^2 + (y4 - y3)^2);
    d4 = sqrt((x1 - x4)^2 + (y1 - y4)^2);

    e1 = (xc - x1)^2 + zc^2;
    e2 = (xc - x2)^2 + zc^2;
    e3 = (xc - x3)^2 + zc^2;
    e4 = (xc - x4)^2 + zc^2;

    h1 = (xc - x1) * (yc - y1);
    h2 = (xc - x2) * (yc - y2);
    h3 = (xc - x3) * (yc - y3);
    h4 = (xc - x4) * (yc - y4);

    % source (s,b) and doublet (q,a) influence in panel coordinates

    if d1 < ep
        s1 = 0;
        q1 = 0;
    else
        f = (y2 - y1) * e1 - (x2 - x1) * h1;
        g = (y2 - y1) * e2 - (x2 - x1) * h2;
        q1 = atan2(zc * (x2 - x1) * (f * r2 - g * r1), zc^2 * (x2 - x1)^2 * r1 * r2 + f * g);
        s1 = ((xc - x1) * (y2 - y1) - (yc - y1) * (x2 - x1)) / d1 * log((r1 + r2 + d1) / (r1 + r2 - d1));
    end

    if d2 < ep
        s2 = 0;
        q2 = 0;
    else
        f = (y3 - y2) * e2 - (x3 - x2) * h2;
        g = (y3 - y2) * e3 - (x3 - x2) * h3;
        q2 = atan2(zc * (x3 - x2) * (f * r3 - g * r2), zc^2 * (x3 - x2)^2 * r2 * r3 + f * g);
        s2 = ((xc - x2) * (y3 - y2) - (yc - y2) * (x3 - x2)) / d2 * log((r2 + r3 + d2) / (r2 + r3 - d2));
    end

    if d3 < ep
        s3 = 0;
        q3 = 0;
    else
        f = (y4 - y3) * e3 - (x4 - x3) * h3;
        g = (y4 - y3) * e4 - (x4 - x3) * h4;
        q3 = atan2(zc * (x4 - x3) * (f * r4 - g * r3), zc^2 * (x4 - x3)^2 * r3 * r4 + f * g);
        s3 = ((xc - x3) * (y4 - y3) - (yc - y3) * (x4 - x3)) / d3 * log((r3 + r4 + d3) / (r3 + r4 - d3));
    end

    if d4 < ep
        s4 = 0;
        q4 = 0;
    else
        f = (y1 - y4) * e4 - (x1 - x4) * h4;
        g = (y1 - y4) * e1 - (x1 - x4) * h1;
        q4 = atan2(zc * (x1 - x4) * (f * r1 - g * r4), zc^2 * (x1 - x4)^2 * r4 * r1 + f * g);
        s4 = ((xc - x4) * (y1 - y4) - (yc - y4) * (x1 - x4)) / d4 * log((r4 + r1 + d4) / (r4 + r1 - d4));
    end

    %% add contribution from the 4 sides

    A = -(q1 + q2 + q3 + q4) / 4 / pi; % times doublet strength

    if abs(zc < ep)
        A = 0;
    end

    B = -(s1 + s2 + s3 + s4) / 4 / pi - zc * A; % times source strength
end
