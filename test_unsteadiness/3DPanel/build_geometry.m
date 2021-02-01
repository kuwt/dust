% Update geometry appending new airfoils to ee,rr,ee_te arrays and incrementing nelems,npoints

function [ ee , rr , ee_te , elems , nelems , npoints ] = ...
            build_geometry( airfoil )


% === Read NACA airfoil and build the geometry ===
if ( strcmp(airfoil.airfoil_str(1:4),'NACA') == 0 )
  % Check that the first 4 digits are 'NACA'
  error([' Only NACA airfoils are implemented. ', ...
          'The airfoil string must beging with NACA. STOP.' ])
end

if (     length(airfoil.airfoil_str) == 8 ) % NACA 4-digit airfoil
  M = str2num(airfoil.airfoil_str(5)  ) ;  
  P = str2num(airfoil.airfoil_str(6)  ) ;
  SS= str2num(airfoil.airfoil_str(7:8)) ;
  mirror = 0 ;
  [ x1 , y1 ]  = setAirfoil4( M , P , SS , ...
         airfoil.chord       , airfoil.nChordPanels , ...
         airfoil.refPoint(1) , airfoil.refPoint(2)  , ...
         airfoil.xcRefPoint  , airfoil.theta , mirror ) ;
                                 
elseif ( length(airfoil.airfoil_str) == 9 ) % NACA 5-digit airfoil
  error([' NACA 5-digit airfoils not implemented yet ' ] ) 
else
  error([' Only NACA 4 and 5-digit airfoils are implemented. ', ...
          'airfoil(',ia,').airfoil_str must be 8 or 9-character string'])
end

% === Append new points and elements to node and connectivity arrays ===
rr = [ x1 ; y1 ] ;
ee = [ (1:length(x1)-1) ; (2:length(x1)) ] ;
ee_te = [ 1 , 2*airfoil.nChordPanels ] ;

nelems  = size(ee,2);
npoints = size(rr,2); 

% === Add new elems to the elems structure ===
% elems is an array of "elements type", objects containing all the relevant info
for ie = 1 : nelems
   elems(ie).airfoilId = airfoil.id ;
   elems(ie).id   = ie ;
   elems(ie).ver1 = rr(:,ee(1,ie)) ;
   elems(ie).ver2 = rr(:,ee(2,ie)) ;
   elems(ie).cen  = 0.5 * ( elems(ie).ver1 + elems(ie).ver2 ) ;
   elems(ie).len  =   norm( elems(ie).ver2 - elems(ie).ver1 ) ;
   elems(ie).tver = (  elems(ie).ver2 - elems(ie).ver1 ) / elems(ie).len ;
   elems(ie).nver = [ -elems(ie).tver(2) ; elems(ie).tver(1) ] ;
end


