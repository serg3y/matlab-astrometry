function [Long,Lat]=pr_ignomonic(X,Y,CenterVec);
%------------------------------------------------------------------------------
% pr_ignomonic function                                               AstroMap
% Description: Project coordinates (X and Y) using the
%              inverse Gnomonic non conformal projection,
% Input  : - Vector of X, in radians.
%          - Vector of Y, in radians.
%          - Central coordinate vector [Long_center,Lat_center],
% Output : - Vector of longitude in radians.
%          - Vector of latitude in radians.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                      June 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: pr_gnomonic.m
% Reliable: 2
%------------------------------------------------------------------------------

Rho   = sqrt(X.^2+Y.^2);
C     = atan(Rho);

Long1 = CenterVec(1);
Lat1  = CenterVec(2);

if (Rho==0),
   Lat = Lat1;
   Long = Long1;
else
   Lat   = asin(cos(C).*sin(Lat1) + Y.*sin(C).*cos(Lat1)./Rho);
   Long  = Long1 + atan(X.*sin(C)./(Rho.*cos(Lat1).*cos(C) - Y.*sin(Lat1).*sin(C)));
end
