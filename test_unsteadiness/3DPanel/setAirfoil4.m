function [x,y] = setAirfoil4(M,P,SS,c,n, x0,y0,ca,theta,mirror )

% percentuale di corda attorno alla quale viene ruotato il profilo
[a,b] = naca4digit(M,P,SS,c,n);

if ( mirror == 1 )
  b = -b;
  a(1:end) = a(end:-1:1);
  b(1:end) = b(end:-1:1);
end


c1 = ca * c;

x = x0 + c1 + cos(theta).*(a-c1) + sin(theta).*b;
y = y0 - sin(theta).*(a-c1) + cos(theta).*b;
