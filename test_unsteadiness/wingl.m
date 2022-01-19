function [u, v, w] = wingl(x, y, z, ib, jb, qf, gamma)
    % calculate induced velocity at a point (x, y, z) due to longitudinal vorticity
    % distribution gamax(i,j) only (semi-span), in a wing fixed coodinate system + (t.e. unsteady vortex)

    u = 0;
    v = 0;
    w = 0;

    for i = 1:ib %7

        for j = 1:jb %7
            [u2, v2, w2] = vortex(x, y, z, ...
                qf(i, j + 1, 1), qf(i, j + 1, 2), qf(i, j + 1, 3), ...
                qf(i + 1, j + 1, 1), qf(i + 1, j + 1, 2), qf(i + 1, j + 1, 3), ...
                gamma(i, j));
            [u4, v4, w4] = vortex(x, y, z, ...
                qf(i + 1, j, 1), qf(i + 1, j, 2), qf(i + 1, j, 3), ...
                qf(i, j, 1), qf(i, j, 2), qf(i, j, 3), ...
                gamma(i, j));
            u = u + u2 + u4;
            v = v + v2 + v4;
            w = w + w2 + w4;
        end %7

    end %7

    % add influence of latest unsteady wake element

    i = ib; 

    for j = 1:jb %8
        [u3, v3, w3] = vortex(x, y, z, ...
            qf(i + 1, j + 1, 1), qf(i + 1, j + 1, 2), qf(i + 1, j + 1, 3), ...
            qf(i + 1, j, 1), qf(i + 1, j, 2), qf(i + 1, j, 3), ...
            gamma(i, j));
        u = u + u3;
        v = v + v3;
        w = w + w3;
    end %8

end
