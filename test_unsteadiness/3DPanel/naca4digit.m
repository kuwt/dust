function [x,y] = naca4digit(M,P,SS,c,n)

m  = M  / 100;
p  = P  / 10;
t  = SS / 100;

if ( m == 0 ) % p must be .ne. 0 !!!
  p = 1 ;
end

%> === Chord discretization ===
xv = linspace(0.0,c,n+1);            % uniform
% xv = c/2.0 .*(1.0-cos(pi.*xv./c));   % cosine
xv = c .*(1.0-cos(0.5*pi.*xv./c));   % half-cosine

%> === Thickness ===
ytfcn = @(x) 5.*t.*c.*(0.2969.*(x./c).^0.5 - 0.1260.*(x./c) ...
    - 0.3516.*(x./c).^2 + 0.2843.*(x./c).^3 - 0.1015.*(x./c).^4);
% 1015: open TE
% 1036: closed TE
yt = ytfcn(xv);

%> === Mean line ===
yc = zeros(size(xv));

for ii = 1 : n+1
    if xv(ii) <= p*c
        yc(ii) = c*(m/p^2 *(xv(ii)/c) * (2*p - (xv(ii)/c)));
    else
        yc(ii) = c*(m/(1-p)^2 * (1 + (2*p - (xv(ii)/c))*(xv(ii)/c) -2*p));
    end
end

%> === Mean line slope ===
dyc = zeros(size(xv));

for ii = 1 : n+1
    if xv(ii) <= p*c
        dyc(ii) = m/p^2 * 2*(p-xv(ii)/c);
    else
        dyc(ii) = m/(1-p)^2 * 2*(p-xv(ii)/c);
    end
end

%> === Upper and lower sides ===
th = atan2(dyc,1);
xU = xv - yt.*sin(th);
yU = yc + yt.*cos(th);
xL = xv + yt.*sin(th);
yL = yc - yt.*cos(th);

%> Sort points
x = zeros(1,2*n+1);
y = zeros(1,2*n+1);
for ii = 1 : n
   x(ii) = xL(n+2-ii);
   y(ii) = yL(n+2-ii);
end

x(n+1:2*n+1) = xU;
y(n+1:2*n+1) = yU;










