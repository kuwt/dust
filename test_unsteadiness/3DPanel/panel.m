function [c1, c2, c3, t1, t2, t3, v1, v2, v3, s] = panel(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

    % x,y,z panel corner points
    % c,t,v chord-wise, tangential, normal vectors
    % s: source strength

    %% calculate chorwise vector

    a1 = ((x2 + x4) - (x1 +x3)) / 2;
    a2 = ((y2 + y4) - (y1 +y3)) / 2;
    a3 = ((z2 + z4) - (z1 +z3)) / 2;
    aa = sqrt(a1^2 + a2^2 + a3^2);
    c1 = a1 / aa;
    c2 = a2 / aa;
    c3 = a3 / aa;

    %% another vector in this plane
    b1 = x4 - x1;
    b2 = y4 - y1;
    b3 = z4 - z1;
    % normal vector
    v1 = c2 * b3 - c3 * b2;
    v2 = b1 * c3 - c1 * b3;
    v3 = c1 * b2 - c2 * b1;
    vv = sqrt(v1^2 + v2^2 +v3^2);
    v1 = v1 / vv;
    v2 = v2 / vv;
    v3 = v3 / vv;
    % tangential vector
    t1 = v2 * c3 - v3 * c2;
    t2 = c1 * v3 - v1 * c3;
    t3 = v1 * c2 - v2 * c1;
    %calculation of panel area
    e1 = x3 - x1;
    e2 = y3 - y1;
    e3 = z3 - z1;
    f1 = x2 - x1;
    f2 = y2 - y1;
    f3 = z2 - z1;
    %normal areas
    s11 = f2 * b3 - f3 * b2;
    s12 = b1 * f3 - f1 * b3;
    s13 = f1 * b2 - f2 * b1;
    s21 = b2 * e3 - b3 * e2;
    s22 = e1 * b3 - b1 * e3;
    s23 = b1 * e2 - b2 * e1;
    s = 0.5 * (sqrt(s11^2 + s12^2 + s13^2) + sqrt(s21^2 + s22^2 + s23^2));
end
