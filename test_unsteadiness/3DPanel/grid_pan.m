function [qf, qc, cr, ds, sigma, s] = grid_pan(rr, jb, b, xtip, ztip, croot, ctip, dxw, ut, wt)

    ib1 = size(rr, 2);

    for i = 1:ib1
        qf(i, 1, 1) = rr(1, i);
        qf(i, 1, 3) = rr(2, i);
    end

    ib = ib1 - 1;
    ib2 = ib1 + 1;
    jb1 = jb + 1; % jb: number of span wise panels

    % calculate panels corner points: qf(i,j,(x,y,z))

    for j = 1:jb1
        y = b / 2 / jb * (j - 1);
        dxle = xtip * 2 * y / b; % b full span, dxle: local sweep
        dzle = ztip * 2 * y / b;
        chord = croot - (croot - ctip) * 2 * y / b;

        for i = 1:ib1
            qf(i, j, 1) = qf(i, 1, 1) * chord + dxle;
            qf(i, j, 2) = y;
            qf(i, j, 3) = qf(i, 1, 3) * chord + dzle;
        end

        % wake far field points

        qf(ib2, j, 1) = qf(ib1, j, 1) + dxw;% check for unsteady
        qf(ib2, j, 2) = qf(ib1, j, 2);
        qf(ib2, j, 3) = qf(ib1, j, 3);
    end

    %% wing collocation points

    for j = 1:jb

        for i = 1:ib1
            qc(i, j, 1) = (qf(i, j, 1) + qf(i, j + 1, 1) + qf(i + 1, j + 1, 1) + qf(i + 1, j, 1)) / 4;
            qc(i, j, 2) = (qf(i, j, 2) + qf(i, j + 1, 2) + qf(i + 1, j + 1, 2) + qf(i + 1, j, 2)) / 4;
            qc(i, j, 3) = (qf(i, j, 3) + qf(i, j + 1, 3) + qf(i + 1, j + 1, 3) + qf(i + 1, j, 3)) / 4;

            %% computation of chord wise vectors, tangential and normal vectors, panel area and source strength
            [ds(i, j, 1), ds(i, j, 2), ds(i, j, 3), ... %chord-wise
                    ds(i, j, 4), ds(i, j, 5), ds(i, j, 6), ... %tangential
                    ds(i, j, 7), ds(i, j, 8), ds(i, j, 9), ... %normal, area ds(i, j, 10)
                    ds(i, j, 10)] = panel(qf(i, j, 1), qf(i, j, 2), qf(i, j, 3), ... %1
                qf(i + 1, j, 1), qf(i + 1, j, 2), qf(i + 1, j, 3), ... %4
                qf(i, j + 1, 1), qf(i, j + 1, 2), qf(i, j + 1, 3), ... %3
                qf(i + 1, j + 1, 1), qf(i + 1, j + 1, 2), qf(i + 1, j + 1, 3)); %2
            sigma(i, j) = ds(i, j, 7) * ut + ds(i, j, 9) * wt;
        end

    end

    s = 0.5 * b * (croot + ctip);
    %c = c / b;
    ar = b * b / s;

    % transform the 4 panel corner points into panel frame of ref: this is needed later to calculate the influence coefficients

    for j = 1:jb

        for i = 1:ib1
            [cr(i, j, 1), cr(i, j, 2), cr(i, j, 3)] = convert(qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), ...
                qf(i, j, 1), qf(i, j, 2), qf(i, j, 3), ...
                ds(i, j, 1), ds(i, j, 2), ds(i, j, 3), ...
                ds(i, j, 4), ds(i, j, 5), ds(i, j, 6), ...
                ds(i, j, 7), ds(i, j, 8), ds(i, j, 9));

            [cr(i, j, 4), cr(i, j, 5), cr(i, j, 6)] = convert(qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), ...
                qf(i +1, j, 1), qf(i + 1, j, 2), qf(i + 1, j, 3), ...
                ds(i, j, 1), ds(i, j, 2), ds(i, j, 3), ...
                ds(i, j, 4), ds(i, j, 5), ds(i, j, 6), ...
                ds(i, j, 7), ds(i, j, 8), ds(i, j, 9));

            [cr(i, j, 7), cr(i, j, 8), cr(i, j, 9)] = convert(qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), ...
                qf(i + 1, j + 1, 1), qf(i + 1, j + 1, 2), qf(i + 1, j + 1, 3), ...
                ds(i, j, 1), ds(i, j, 2), ds(i, j, 3), ...
                ds(i, j, 4), ds(i, j, 5), ds(i, j, 6), ...
                ds(i, j, 7), ds(i, j, 8), ds(i, j, 9));

            [cr(i, j, 10), cr(i, j, 11), cr(i, j, 12)] = convert(qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), ...
                qf(i, j + 1, 1), qf(i, j + 1, 2), qf(i, j + 1, 3), ...
                ds(i, j, 1), ds(i, j, 2), ds(i, j, 3), ...
                ds(i, j, 4), ds(i, j, 5), ds(i, j, 6), ...
                ds(i, j, 7), ds(i, j, 8), ds(i, j, 9));
        end

    end

end
