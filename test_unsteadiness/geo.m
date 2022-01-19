function [s, qf, qc, ds, ar] = geo(b, c, alf, ib, jb, alamda, bb, dxw, alfa)

    % ib: chordwise boxes
    % jb: spanwise boxes

    ib1 = ib + 1;
    jb1 = jb + 1;

    for i = 1:ib1 %2 
        sn(i) = sin(alf(i));
        cs(i) = cos(alf(i));
    end %2 

    ctg1 = tan(pi / 2 - alamda(1));
    ctg2 = tan(pi / 2 - alamda(2));

    ctip = c + b * (ctg2 - ctg1);
    s = b * (c + ctip) / 2;
    ar = 2 * b^2 / s; % aspect ratio

    %% wing fixed vortices location  (qf(i,j,(x,y,z)...)

    bj = 0; 

    for j = 1:jb1 %3 

        if j > 1 
            bj = bj + bb(j - 1);
        end 

        z1 = 0; 
        dc1 = bj * ctg1; % leading edge x
        dc2 = bj * ctg2; % trailing edge x
        dx1 = (c + dc2 - dc1) / ib;

        for i = 1:ib %1
            qf(i, j, 1) = dc1 + dx1 * (i - 0.75);
            qf(i, j, 2) = bj;
            qf(i, j, 3) = z1 - 0.25 * dx1 * sn(i);
            z1 = z1 - dx1 * sn(i);
        end %1

        % the following lines are due to wake distance from trailing edge

        qf(ib1, j, 1) = c + dc2 + dxw;
        qf(ib1, j, 2) = qf(ib, j, 2);
        qf(ib1, j, 3) = z1 - dxw * sn(ib);

    end %3 

    %% wing collocation points

    for j = 1:jb %4
        z1 = 0;
        bj = qf(1, j, 2) + bb(j) / 2;
        dc1 = bj * ctg1;
        dc2 = bj * ctg2;
        dx1 = (c + dc2 - dc1) / ib;

        for i = 1:ib %4
            qc(i, j, 1) = dc1 + dx1 * (i - 0.25);
            qc(i, j, 2) = bj;
            qc(i, j, 3) = z1 - 0.75 * dx1 * sn(i);
            z1 = z1 - dx1 * sn(i);
            ds(i, j) = dx1 * bb(j);
        end %4

    end %4

    %% rotation of wing points due to alpha

    sn1 = sin(-alfa);
    cs1 = cos(-alfa);

    for i = 1:ib1 %6

        for j = 1:jb1 %6

            qf1 = qf(i, j, 1); 
            qf(i, j, 1) = qf1 * cs1 - qf(i, j, 3) * sn1;
            qf(i, j, 3) = qf1 * sn1 + qf(i, j, 3) * cs1;

            if i == ib1 || j >= jb1 
            else

                qc1 = qc(i, j, 1);
                qc(i, j, 1) = qc1 * cs1 - qc(i, j, 3) * sn1;
                qc(i, j, 3) = qc1 * sn1 + qc(i, j, 3) * cs1;

            end
            % continue 6
        end %6

    end %6 

end
