% 3D-Panel code for simple wing planforms
clc; clear all; close all;
%% Input data

alpha1 = 5.0; % AoA
croot = 1; % root chord
ctip = 1; % tip chord
xtip = 0; % aft sweep of tip
ztip = 0; % dihedro
b = 10; % wing span
vt = 10; % free stream speed

jb = 13; % span-wise counters
dxw = 100 * b; % wake (fix for unsteadiness)
ro = 1.0;

ut = vt * cos(alpha1 * pi / 180);
wt = vt * sin(alpha1 * pi / 180);

airfoil.id = 1;
airfoil.airfoil_str = 'NACA0012';
airfoil.chord = 1.0;
airfoil.theta = deg2rad(0.0);
airfoil.xcRefPoint = 0.0;
airfoil.refPoint = [0; 0.0];
airfoil.nChordPanels = 19;
ib = airfoil.nChordPanels; % points
%% wing geometry
[ee, rr, ee_te, elems, nelems, npoints] = build_geometry(airfoil);
% close trailing edge
rr(2, 1) = 0;
rr(2, end) = 0;

[qf, qc, cr, ds, sigma, s] = grid_pan(rr, jb, b, xtip, ztip, croot, ctip, dxw, ut, wt);

mesh(qf(:, :, 1), qf(:, :, 2), qf(:, :, 3), 'edgecolor', 'k')
axis equal

ib1 = ib + 1;
%ib2 = ib + 2;
jb1 = jb + 1;

%%  aerodynamic calculation

% influence coefficients calculation
% collocation point counter

k = 0;

for i = 1:ib

    for j = 1:jb
        k = k + 1;
        l = 0;
        rh = 0;
        % influencing panel counter
        for i1 = 1:ib

            for j1 = 1:jb
                l = l + 1;

                if i1 == 1
                    % calculate wake contribution: first convert collocation point to panel coordinates
                    % and then calculate influence coefficients
                    [xc, yc, zc] = convert(qc(ib1, j1, 1), qc(ib1, j1, 2), qc(ib1, j1, 3), ...
                        qc(ib1, j1, 1), qc(ib1, j1, 2), qc(ib1, j1, 3), ...
                        ds(ib1, j1, 1), ds(ib1, j1, 2), ds(ib1, j1, 3), ...
                        ds(ib1, j1, 4), ds(ib1, j1, 5), ds(ib1, j1, 6), ...
                        ds(ib1, j1, 7), ds(ib1, j1, 8), ds(ib1, j1, 9));
                    [wdub, dsig] = influence(xc, yc, zc, ...
                    cr(ib1, j1, 1), cr(ib1, j1, 2), cr(ib1, j1, 3), ...
                        cr(ib1, j1, 4), cr(ib1, j1, 5), cr(ib1, j1, 6), ...
                        cr(ib1, j1, 7), cr(ib1, j1, 8), cr(ib1, j1, 9), ...
                        cr(ib1, j1, 10), cr(ib1, j1, 11), cr(ib1, j1, 12));

                    % add wing's image (symmetry is assumed)
                    [xc, yc, zc] = convert(qc(ib1, j1, 1), qc(ib1, j1, 2), qc(ib1, j1, 3), ...
                        qc(i, j, 1), -qc(i, j, 2), qc(i, j, 3), ...
                        ds(ib1, j1, 1), ds(ib1, j1, 2), ds(ib1, j1, 3), ...
                        ds(ib1, j1, 4), ds(ib1, j1, 5), ds(ib1, j1, 6), ...
                        ds(ib1, j1, 7), ds(ib1, j1, 8), ds(ib1, j1, 9));
                    [wdub1, dsig1] = influence(xc, yc, zc, ...
                    cr(ib1, j1, 1), cr(ib1, j1, 2), cr(ib1, j1, 3), ...
                        cr(ib1, j1, 4), cr(ib1, j1, 5), cr(ib1, j1, 6), ...
                        cr(ib1, j1, 7), cr(ib1, j1, 8), cr(ib1, j1, 9), ...
                        cr(ib1, j1, 10), cr(ib1, j1, 11), cr(ib1, j1, 12));

                    ddubj(j1) = wdub + wdub1;
                    dmu2 = ddubj(j);

                else

                    dmu2 = ddubj(j1);
                end

                if i1 == ib
                    dmu2 = -ddubj(j1);
                end % end of wake influence calculation

                % convert collocation point to panel coordinates

                [xc, yc, zc] = convert(qc(i1, j1, 1), qc(i1, j1, 2), qc(ib1, j1, 3), ...
                    qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), ...
                    ds(i1, j1, 1), ds(i1, j1, 2), ds(i1, j1, 3), ...
                    ds(i1, j1, 4), ds(i1, j1, 5), ds(i1, j1, 6), ...
                    ds(i1, j1, 7), ds(i1, j1, 8), ds(i1, j1, 9));

                [dmu, dsig] = influence(xc, yc, zc, ... 
                cr(i1, j1, 1), cr(i1, j1, 2), cr(i1, j1, 3), ...
                    cr(i1, j1, 4), cr(i1, j1, 5), cr(i1, j1, 6), ...
                    cr(i1, j1, 7), cr(i1, j1, 8), cr(i1, j1, 9), ...
                    cr(i1, j1, 10), cr(i1, j1, 11), cr(i1, j1, 12));

                if i1 == i && j1 == j
                    dmu = -0.5; % a panel influence on itself is dmu = 1/2
                end

                % add influence of wing's image (other half)

                [xc, yc, zc] = convert(qc(i1, j1, 1), qc(i1, j1, 2), qc(ib1, j1, 3), ...
                    qc(i, j, 1), -qc(i, j, 2), qc(i, j, 3), ...
                    ds(i1, j1, 1), ds(i1, j1, 2), ds(i1, j1, 3), ...
                    ds(i1, j1, 4), ds(i1, j1, 5), ds(i1, j1, 6), ...
                    ds(i1, j1, 7), ds(i1, j1, 8), ds(i1, j1, 9));

                [dmu1, dsig1] = influence(xc, yc, zc, ...
                cr(i1, j1, 1), cr(i1, j1, 2), cr(i1, j1, 3), ...
                    cr(i1, j1, 4), cr(i1, j1, 5), cr(i1, j1, 6), ...
                    cr(i1, j1, 7), cr(i1, j1, 8), cr(i1, j1, 9), ...
                    cr(i1, j1, 10), cr(i1, j1, 11), cr(i1, j1, 12));

                A(k, l) = dmu + dmu1 - dmu2; 

                rh = rh + (dsig + dsig1) * sigma(i1, j1);
            end

        end

        % calculate rhs

        rhs(k) = rh;

    end

end

k1 = ib * jb;

for k = 1:k1
    dub1(k) = rhs(k);
end

A_inv = A \ eye(size(A, 1)); % non variable wing geometry with time: inversion performed only once
dub1 = A_inv * dub1'; % 13.147
dub1 = dub1';

% wing doublet lattice listing

k = 1;

for i = 1:ib

    for j = 1:jb
        dub(i, j) = dub1(k); 
    end

end

for j = 1:jb
    dub(ib1, j) = dub(1, j) - dub(ib, j);
end

%% force calculation

fl = 0;
fd = 0;
fm = 0;

que = 0.5 * ro * vt^2;

for j = 1:jb

    for i = 1:ib
        i1 = i - 1;
        i2 = i + 1;
        j1 = j - 1;
        j2 = j + 1;

        if i == 1
            i1 = 1;
        end

        if i == ib
            i2 = ib;
        end

        if j == 1
            j1 = 1;
        end

        if j == jb
            j2 = jb;
        end

        % chord wise velocity

        xf = 0.5 * (qf(i + 1, j, 1) + qf(i + 1, j + 1, 1));
        yf = 0.5 * (qf(i + 1, j, 2) + qf(i + 1, j + 1, 2));
        zf = 0.5 * (qf(i + 1, j, 3) + qf(i + 1, j + 1, 3));
        xr = 0.5 * (qf(i, j, 1) + qf(i, j + 1, 1));
        yr = 0.5 * (qf(i, j, 2) + qf(i, j + 1, 2));
        zr = 0.5 * (qf(i, j, 3) + qf(i, j + 1, 3));

        dx2 = qc(i2, j, 1) - xf;
        dy2 = qc(i2, j, 2) - yf;
        dz2 = qc(i2, j, 3) - zf;
        dx3 = qc(i1, j, 1) - xr;
        dy3 = qc(i1, j, 2) - yr;
        dz3 = qc(i1, j, 3) - zr;

        dl1 = sqrt((xf - xr)^2 + (yf - yr)^2 + (zf - zr)^2);
        dl2 = sqrt(dx2^2 + dy2^2 + dz2^2);
        dl3 = sqrt(dx3^2 + dy3^2 + dz3^2);
        dll = dl1 + dl2 + dl3;

        if i == 1
            dll = dl1 / 2 + dl2;
        end

        if i == ib
            dll = dl1 / 2 + dl3;
        end

        ql = -(dub(i2, j) - dub(i1, j)) / dll;

        % spanwise velocity

        dx = qc(i, j2, 1) - qc(i, j1, 1);
        dy = qc(i, j2, 2) - qc(i, j1, 2);
        dz = qc(i, j2, 3) - qc(i, j1, 3);
        dr = sqrt(dx^2 + dy^2 + dz^2);

        qm = -(dub(i, j2) - dub(i, j1)) / dr;

        %first order correction for panel sweep

        ql = ql + qm * (dx^2 + dz^2) / dr;
        qm = qm * (dy^2 + dz^2) / dr;
        qinf = ut * ds(i, j, 9) - wt * ds(i, j, 7);
        cp(i, j) = 1 - ((qinf + ql)^2 + qm^2) / vt^2;
        dl(i, j) = -cp(i, j) * ds(i, j, 10) * ds(i, j, 9);
        dd(i, j) = cp(i, j) * ds(i, j, 10) * ds(i, j, 7);

        fl = fl + dl(i, j);
        fd = fd + dd(i, j);
        fm = fm + dl(i, j) * qc(i, j, 1);
    end

end

cl = fl / (que * s);
cd = fd / (que * s);
cm = fm / (que * s * croot);
