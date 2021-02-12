function [u, v, w] = veloce(x, y, z, it, js1, js2, cs1, sn1, qw, vortic, jb, gamma, qf, ib, sign, sn0, cs0, ch, sx, sz)

    x1 = (x - sx) * cs1 + (z - sz) * sn1;
    y1 = y;
    z1 = -(x - sx) * sn1 + (z - sz) * cs1;

    [u1, v1, w1] = wake(x, y, z, vortic, qw, it, jb);
    [u2, v2, w2] = wake(x, -y, z, vortic, qw, it, jb);

    [u3, v3, w3] = wing(x1, y1, z1, gamma, qf, ib, jb, sign, sn0, cs0);
    [u4, v4, w4] = wing(x1, -y1, z1, gamma, qf, ib, jb, sign, sn0, cs0);

    u33 = cs1 * (u3 + u4) - sn1 * (w3 + w4);
    w33 = sn1 * (u3 + u4) + cs1 * (w3 + w4);

    % influence of mirror image

    if ch > 100
        u5 = 0;
        u6 = 0;
        u77 = 0;
        v5 = 0;
        v6 = 0;
        v7 = 0;
        v8 = 0;
        w5 = 0;
        w6 = 0;
        w77 = 0;
    else % ground effect
        x2 = (x - sx) * cs1 + (-z -sz) * sn1;
        z2 = -(x - sx) * sn1 + (-z -sz) * cs1;
        [u5, v5, w5] = wake (x, y, -z, vortic, qw, it, jb);
        [u6, v6, w6] = wake (x, -y, -z, vortic, qw, it, jb);

        [u7, v7, w7] = wing(x2, y1, z2, gamma, qf, ib, jb, sign, sn0, cs0);
        [u8, v8, w8] = wing(x2, -y1, z2, gamma, qf, ib, jb, sign, sn0, cs0);

        u77 = cs1 * (u7 + u8) - sn1 * (w7 + w8);
        w77 = sn1 * (u7 + u8) + cs1 * (w7 + w8);
    end

    u = u1 + u2 + u33 + u5 + u6 + u77;
    v = v1 - v2 + v3 - v4 + v5 - v6 + v7 - v8;
    w = w1 + w2 + w33 - w5 - w6 - w77;

end
