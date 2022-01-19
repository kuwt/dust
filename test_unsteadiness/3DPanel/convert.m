function [xp,yp,zp] = convert(x0,y0,z0,xb,yb,zb,c1,c2,c3,t1,t2,t3,v1,v2,v3)

    % transformation of a field point b into panel coordinates 0 

    xp = (xb - x0) * c1 + (yb - y0) * c2 + (zb - z0) * c3;
    yp = (xb - x0) * t1 + (yb - y0) * t2 + (zb - z0) * t3;
    zp = (xb - x0) * v1 + (yb - y0) * v2 + (zb - z0) * v3;

end