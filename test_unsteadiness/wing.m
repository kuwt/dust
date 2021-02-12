% CALCULATES INDUCED VELOCITY AT A POINT (X,Y,Z), DUE TO VORTICITY
% DISTRIBUTION GAMA(I,J), OF SEMI-CONFIGURATION - IN A WING FIXED
% COORDINATE SYSTEM.
function [u, v, w, A1] = wing(x, y, z, gamma, qf, ib, jb, sign, sn0, cs0)

    u = 0;
    v = 0;
    w = 0;

    for i = 1:ib %7 

        for j = 1:jb %7
            [u1, v1, w1] = vortex(x, y, z, ...
                qf(i, j, 1), qf(i, j, 2), qf(i, j, 3), ...
                qf(i, j + 1, 1), qf(i, j + 1, 2), qf(i, j + 1, 3), ...
                gamma(i, j));
            [u2, v2, w2] = vortex(x, y, z, ...
                qf(i, j + 1, 1), qf(i, j + 1, 2), qf(i, j + 1, 3), ...
                qf(i + 1, j + 1, 1), qf(i + 1, j + 1, 2), qf(i + 1, j + 1, 3), ...
                gamma(i, j));
            [u3, v3, w3] = vortex(x, y, z, ...
                qf(i + 1, j + 1, 1), qf(i + 1, j + 1, 2), qf(i + 1, j + 1, 3), ...
                qf(i + 1, j, 1), qf(i + 1, j, 2), qf(i + 1, j, 3), ...
                gamma(i, j));
            [u4, v4, w4] = vortex(x, y, z, ...
                qf(i + 1, j, 1), qf(i + 1, j, 2), qf(i + 1, j, 3), ...
                qf(i, j, 1), qf(i, j, 2), qf(i, j, 3), ...
                gamma(i, j));

            u0 = u1 + u2 + u3 + u4;
            v0 = v1 + v2 + v3 + v4;
            w0 = w1 + w2 + w3 + w4;
            A1(i, j) = u0 * sn0(i) + w0 * cs0(i);

            if sign >= 1
                A1(i, j) = u0 * sn0(i) - w0 * cs0(i);
            end

            u = u + u0;
            v = v + v0;
            w = w + w0;
        end %7

    end %7

end
