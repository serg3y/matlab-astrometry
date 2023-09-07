function [Long,Lat] = xy2sky_tan(FileName,X,Y)
%Given a FITS image, SIM or a structure containing the FITS WCS keyword
%(returned by fits_get_wcs.m), where the WCS is represented using the
%tangential projection, convert X and Y position in the image to longitude
%and latitude.
% Input  : - A string containing a FITS image name, SIM, or a structure
%            containing the FITS WCS keyword, returned by fits_get_wcs.m.
%          - Vector of X positions [pixel].
%          - Vector of Y positions [pixel].
%          - HDU number in FITS image from which to read the header.
%            Default is 1.
% Output : - Column vector of Longitudes [radian].
%          - Column vector of Latitude [radian].
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Nov 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%Example: 
% [X,Y] = xy2sky_tan('File.fits',1,1);
% WCS = get_fits_wcs('File.fits');
% [X,Y] = xy2sky_tan(WCS,1,1);
% Sim = images2sim('File.fits');
% [X,Y] = xy2sky_tan(Sim,1,1);

% deal with types of input
if isstruct(FileName)
    if isfield(FileName,'CTYPE1')
        % assume input is a WCS structure (generated by fits_get_wcs.m)
        WCS = FileName;
        CallGetWCS = false;
    else
        CallGetWCS = true;
    end
else
    CallGetWCS = true;
end

if CallGetWCS
    error('Not implemented')
end

% transformation
if ~strcmp(WCS.CUNIT1,WCS.CUNIT2)
     error('CUNIT1 must be identical to CUNIT2')
end

switch lower(WCS.CUNIT1)
    case {'deg' 'degree'}
       Factor = 180/pi;
    case 'rad'
       Factor = 1;
    otherwise
       error('Unknown CUNIT option')
end

if size(X,2)==1
   X = X';
end

if size(Y,2)==1
   Y = Y';
end

D = WCS.CD*[X - WCS.CRPIX1; Y - WCS.CRPIX2]; % + [WCS.CRVAL1; WCS.CRVAL2]
[Long,Lat] = pr_ignomonic(D(1,:)'/Factor,D(2,:)'/Factor,[WCS.CRVAL1 WCS.CRVAL2]/Factor);