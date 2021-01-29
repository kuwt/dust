clear; close all; clc
%
%% unsteady rectangular lifting surface (vlm)
%---------------------------------------------------------
% this is a 3-d linear code for rectangular planforms (with ground effect)
% in unsteady motion using the vortex lattice method (by Joe Katz, 1975)
%
%% MODES OF OPERATION
% 1. steady state % : set dt=dx/vt*ib*10, and nsteps=5
% 2. sudden acceleration : set dt=dx/vt/4. and nsteps= up to 50
% 3. heaving oscillations: bh=heaving ampl. om= freq.
% 4. pitch oscillations : omega=frequency, teta=momentary angle
% 5. for computational economy the parameter iw might be used
% note; induced drag calculation increases computation time and
% can be disconnected.

%% input data
tic 
% geometry 
ib = 4; % number of chordwise panels
jb = 13; % no. of spanwise panels

nw = 2; % nw - is the number of (timewise) deforming wake elements.
ro = 1.225; % density
bh = 0.; % heaving amplitude
om = 0*(2*pi); % frequency
theta = 0; % pitch
omega = 0.;
vt = 50.0; % far field velocity
c = 1.; % chord
b = 2; % semi span
dx = c / ib; % chord panel dimensions
dy = b / jb; % span panel dimensions
ch = 10000 .* c; % ground clearance
alfa1 = 5;
alfa0 = 5;

alfa = (alfa1 + alfa0) * pi / 180.0;
alamda(1) = 90 * pi / 180; % sweep back angles. (alamda < 90, sweep backward).
alamda(2) = alamda(1);
w11 = 0; 

%% time settings 
tf = 0.2; % final time 
%dt=dx/vt/4; % unsteady 
%dt = dx / vt *ib*10; % steady 
dt = 1.25e-3; 
time = 0:dt:tf;
nsteps = numel(time); 


%%  Data initialization
%alf = zeros(size(nw,1));
%sn0 = zeros(size(nw,1));
%cs0 = zeros(size(nw,1));
%alam = zeros(size(jb,1));

for i = 1:ib%2
    alf(i) = 0.;
end %2

alf(ib + 1) = alf(ib);

t = -dt; % time in seconds
dxw = 0.3 * vt * dt;

for j = 1:jb%3
    bb(j) = dy;
end %3

%% constants

k = 0; % fa schifo (meglio ones, zeros)

for i = 1:ib %1

    for j = 1:jb%1
        k = k + 1;
        ww(k) = 0.0;
        dlt(i, j) = 0.0;
        vortic(i, j) = 0.0;
        vort1(i, j) = 0.0;
        gamma(i, j) = 1.; % is required for influence matrix calculations.
    end %1

end %1

%% calculation of collocation points.
%
% geo calculates wing collocation points qc, and vortex tips qf, normal
% ds, surface s
[s, qf, qc, ds, ar] = geo(b, c, alf, ib, jb, alamda, bb, dxw, alfa);

alam1 = alamda(1) * 180 / pi;
alam2 = alamda(2) * 180 / pi;

for i = 1:ib%31
    all = alf(i) * 180 / pi;
end %31

ib1 = ib + 1;
jb1 = jb + 1;

% TODO: add interpolator
figure 
%plot3(qf(:,:,1), qf(:,:,2), qf(:,:,3), '*')
hold on
%plot3(qc(:,:,1), qc(:,:,2), qc(:,:,3), 'o')
grid on 
axis equal
%% PROGRAM START

for it = 1:nsteps%100

    t = t + dt;

    % path information (add aeroelasticity)

    sx = -vt * t; % position x
    dsx = -vt; % velocity x
    ch1 = ch;

    if ch > 100.0
        ch1 = 0;
    end

    sz = bh * sin(om * t) + ch1; % position z  (plunge)
    dsz = bh * om * cos(om * t); % velocity z

    vt = -cos(theta) * dsx - sin(theta) * dsz; % check!! 
    sn1 = sin(theta);
    cs1 = cos(theta);
    wt = sn1 * dsx - cs1 * dsz;

    for i = 1:ib    %6
        sn0(i) = sin(alfa + alf(i));
        cs0(i) = cos(alfa + alf(i));
    end             %6 

    %% vortex wake shedding points

    for j = 1:jb1 %7
        qw(it, j, 1) = qf(ib1, j, 1) * cs1 - qf(ib1, j, 3) * sn1 + sx;
        qw(it, j, 2) = qf(ib1, j, 2);
        qw(it, j, 3) = qf(ib1, j, 1) * sn1 + qf(ib1, j, 3) * cs1 + sz;
    end         %7
    
    %% aerodynamic calculations
    hold on 
    % influence coefficients calculation

    k = 0; 

    for i = 1:ib %14

        for j = 1:jb %14 
            sign = 0;
            k = k + 1;

            if it > 1  % goto to 12 
            else
                % matrix coefficients calculation occurs only once for the time - fixed - geometry wing. (hinge ?)
                [u, v, w, A1] = wing(qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), gamma, qf, ib, jb, sign, sn0, cs0);
                l = 0;

                for i1 = 1:ib %10

                    for j1 = 1:jb %10
                        l = l + 1;
                        % A(k,l) is the normal velocity component due to a unit vortex lattice
                        A(k, l) = A1(i1, j1);
                    end %10

                end %10

                % add influence of wing other half part
                [u, v, w, A1] = wing(qc(i, j, 1), -qc(i, j, 2), qc(i, j, 3), gamma, qf, ib, jb, sign, sn0, cs0);
                l = 0;

                for i1 = 1:ib %11

                    for j1 = 1:jb %11
                        l = l + 1;
                        % A(k,l) is the normal velocity component due to a unit vortex lattice
                        A(k, l) = A(k, l) + A1(i1, j1);
                    end %11

                end %11

                if ch > 100 
                else
                    % add influence of mirror image
                    sign = 10;  
                    xx1 = qc(i, j, 1) * cs1 - qc(i, j, 3) * sn1 + sx;
                    zz1 = qc(i, j, 1) * sn1 + qc(i, j, 3) * cs1 + sz;
                    xx2 = (xx1 - sx) * cs1 + (-zz1 - sz) * sn1;
                    zz2 = -(xx1 - sx) * sn1 + (-zz1 -sz) * cs1;

                    [u, v, w, A1] = wing(xx2, qc(i, j, 2), zz2, gamma, qf, ib, jb, sign, sn0, cs0);

                    l = 0;

                    for i1 = 1:ib %8

                        for j1 = 1:jb %8
                            l = l + 1;
                            % A(k,l) is the normal velocity component due to a unit vortex lattice
                            A(k, l) = A(k, l) + A1(i1, j1);
                        end %8

                    end %8 

                    % add mirror image influence of wing's other half
                    [u, v, w, A1] = wing(xx2, -qc(i, j, 2), zz2, gamma, qf, ib, jb, sign, sn0, cs0); 

                    l = 0;

                    for i1 = 1:ib %9

                        for j1 = 1:jb %9
                            l = l + 1;
                            % A(k,l) is the normal velocity component due to a unit vortex lattice
                            A(k, l) = A(k, l) + A1(i1, j1);
                        end %9

                    end %9

                    sign = 0;
                end

            end 
            % 12 continue 
            if it == 1 
            else
                % calculate wake influence
                xx1 = qc(i, j, 1) * cs1 - qc(i, j, 3) * sn1 + sx;
                zz1 = qc(i, j, 1) * sn1 + qc(i, j, 3) * cs1 + sz;

                [u, v, w] = wake (xx1, qc(i, j, 2), zz1, vortic, qw, it, jb);
                [u1, v1, w1] = wake (xx1, -qc(i, j, 2), zz1, vortic, qw, it, jb);

                if ch > 100 
                    u2 = 0; v2 = 0; w2 = 0; %121
                    u3 = 0; v3 = 0; w3 = 0;
                else
                    [u2, v2, w2] = wake (xx1, qc(i, j, 2), -zz1, vortic, qw, it, jb);
                    [u3, v3, w3] = wake (xx1, -qc(i, j, 2), -zz1, vortic, qw, it, jb);
                end
                % continue 122
                % wake induced velocity is given in inertial frame
                u = u + u1 + u2 + u3;
                w = w + w1 - w2 - w3;
                u11 = u * cs1 + w * sn1;
                w11 = -u * sn1 + w * cs1;

                % ww(k) is the normal component of the wake
                ww(k) = u11 * sn0(i) + w11 * cs0(i);
            end

            %% calculate wing geometrical downwash
            dw(k) = -vt * sn0(i) + qc(i, j, 1) * omega - wt; % for general motion dw(k) = -vt*sin(alfa) + omega * x

            wts(i, j) = w11; % w11 is positive since the latest unsteady wake element is included in subroutine wing

        end %14

    end %14

    %% solution of the problem: dw(i) = ww(i) + A(i,j) * gamma(i)

    k1 = ib * jb; 

    for k = 1:k1 % 15
        gamma1(k) = dw(k) - ww(k); % rhs
    end %15

    if it > 1
    else
        A_inv = A \ eye(size(A, 1)); % non variable wing geometry with time: inversion performed only once
    end

    gamma1 = A_inv * gamma1'; % 13.147
    gamma1 = gamma1';
    % wing vortex lattice listing

    k = 0; 

    for i = 1:ib %17

        for j = 1:jb %17
            k = k + 1;
            gamma(i, j) = gamma1(k);
        end

    end

    %% wake shedding

    for j = 1:jb %171
        vortic(it, j) = gamma(ib, j);
        vortic(it + 1, j) = 0;
    end %171

    %% Wake rollup calculation

    iw = 1;

    if it == 1 % goto 193 
    else
   
        if it >= nw
            iw = it - nw + 1; % nw is the number of timewise deforming wake elemnts (n_dt_wake)
        end

        i1 = it - 1;  
        js1 = 0; % ?
        js2 = 0;

        for i = iw:i1 %18

            for j = 1:jb1 %18 

                [u, v, w] = veloce(qw(i, j, 1), qw(i, j, 2), qw(i, j, 3), it, js1, js2, cs1, sn1, qw, vortic, jb, gamma, qf, ib, sign, sn0, cs0, ch, sx, sz);
                uvw(i, j, 1) = u * dt; 
                uvw(i, j, 2) = v * dt;
                uvw(i, j, 3) = w * dt;
            end %19

        end %19

        for i = iw:i1 %19

            for j = 1:jb1 %19
                qw(i, j, 1) = qw(i, j, 1) + uvw(i, j, 1);
                qw(i, j, 2) = qw(i, j, 2) + uvw(i, j, 2);
                qw(i, j, 3) = qw(i, j, 3) + uvw(i, j, 3);
            end %19

        end %19
%         plot3(qw(iw,:,1), qw(iw,:,2), qw(iw,:,3),'*')
%         xlabel('x')
%         ylabel('y')
%         zlabel('z')
%         drawnow
%         view([45 45]);
    end 
    % continue 193 

    %% Force calculation

    fl = 0;
    fd = 0;
    fm = 0;
    fg = 0;
    que = 0.5 * ro * vt * vt; 

    for j = 1:jb %20 
        sigma = 0;
        sigma1 = 0;
        dly(j) = 0;

        for i = 1:ib %20 

            if i == 1 
                gammaij = gamma(i, j);
            end

            if i > 1
                gammaij = gamma(i, j) - gamma(i - 1, j);
            end

            dxm = (qf(i, j, 1) + qf(i, j + 1, 1)) / 2 ; % distance from the leading edge
            dxm = dxm - 0.25*c; % shift at 0.25
            sigma1 = (0.5 * gammaij + sigma) * dx;

            sigma = gamma(i, j);
            dfdt = (sigma1 - dlt(i, j)) / dt; % velocity potential time derivative

            dlt(i, j) = sigma1;
            dl(i, j) = ro * (vt * gammaij + dfdt) * bb(j) * cs0(i);

            %% induced drag calculation
            [u1, v1, w1] = wingl(qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), ib, jb, qf, gamma);
            [u2, v2, w2] = wingl(qc(i, j, 1), -qc(i, j, 2), qc(i, j, 3), ib, jb, qf, gamma);

            if ch > 100
                w3 = 0; %194
                w4 = 0;
            else
                xx1 = qc(i, j, 1) * cs1 - qc(i, j, 3) * sn1 + sx;
                zz1 = qc(i, j, 1) * sn1 + qc(i, j, 3) * cs1 + sz;
                xx2 = (xx1 -sx) * cs1 + (-zz1 - sz) * sn1;
                zz2 = -(xx1 -sx) * sn1 + (-zz1 - sz) * cs1;
                [u3, v3, w3] = wingl(qc(i, j, 1), qc(i, j, 2), qc(i, j, 3), ib, jb, qf, gamma);
                [u4, v4, w4] = wingl(qc(i, j, 1), -qc(i, j, 2), qc(i, j, 3), ib, jb, qf, gamma);
            end

            w8 = w1 + w2 - w3 - w4; %195 

            % add influence of mirror image (ground)

            cts = -(wts(i, j) + w8) / vt;
            dd1 = ro * bb(j) * dfdt * sn0(i);
            dd2 = ro * bb(j) * vt * gammaij * cts;
            dd(i, j) = dd1 + dd2;

            dp(i, j) = dl(i, j) / ds(i, j) / que;
            dly(j) = dly(j) + dl(i, j);
            fl = fl + dl(i, j);
            fd = fd + dd(i, j);
            fm = fm + dl(i, j) * dxm;
            fg = fg + gammaij*bb(j); 
        end  %20

    end %20

    cl(it) = fl / (que * s*2);
    cd = fd / (que * s*2);
    cm(it) = fm / (que * s*2 * c);
    cl00 = 2 * pi * alfa / (1 + 2 / ar);

    if abs(cl00) < 1e-20
        cl00 = cl(it);
    end

    clt = cl(it) / cl00;
    cfg = fg / (0.5 * vt * s) / cl00;
    
end %100

%%
figure  
plot(time*vt/c, cl)
title('CL')
xlabel('U t/c')
ylim([0,1])
toc 


