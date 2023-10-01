function setup_Astrometry(fov)
%Auto install WSL, Ubuntu, Astrometry.net for Windows 10.
% setup_Astrometry(fov)
%
%Manual install:
%-Install Ubuntu using windows command line:
% wsl --install           %installs the default distro (probably Ubuntu)
% wsl --install Ubuntu -n %install Ubuntu but do not launch it
% Note: If windows does ont have update 22H2 it may not have WSL feature.
% Note: "wsl --install" may need to be run before other commands work.
%-Reboot computer after installing Ubuntu
%-Install astrometry.net, from Ubuntu run:
% sudo apt udate                    %required for next step
% sudo apt install astrometry.net   %install astrometry software
%-Download index files down to 10% of the FOV, eg 5' for 0.8° FOV:
% sudo apt install astrometry.net astrometry-data-2mass-08-19  %156MB
% sudo apt install astrometry.net astrometry-data-2mass-07     %161MB
% sudo apt install astrometry.net astrometry-data-2mass-06     %328MB
% sudo apt install astrometry.net astrometry-data-2mass-05     %659MB
% Note: This downloads the 4200-series from: http://data.astrometry.net
%       The 5200-Lite may be better, its based on GAIA DR2 and Tycho2
%       The 5200-Heave includes G/BP/RP mags, proper motions, parallaxes
%-Install Source-Extractor software (optional, not used at the moment):
% sudo apt install sextractor
%-Test (optional):
% !bash -c "solve-field /mnt/c/MatLab/matlabtoolbox.git/data/DSTG_VIS_500/32711_NAVSTAR_62_USA_201/20201214_143228.fit --overwrite --downsample 2"
% !wsl      solve-field /mnt/c/MatLab/matlabtoolbox.git/data/DSTG_VIS_500/32711_NAVSTAR_62_USA_201/20201214_143228.fit --overwrite --downsample 2
%-Ref:
% https://www.hnsky.org/linux_subsyst.htm
% https://github.com/Jusas/astrometry-api-lite#installation
% https://learn.microsoft.com/en-us/windows/wsl/install
% https://packages.debian.org/source/bookworm/astrometry-data-2mass
% http://data.astrometry.net/

% Defaults
if nargin<2 || ismept(fov), fov = 0.86; end
fold = fileparts(mfilename('fullpath'));

% Install WSL
fprintf('Installing WSL\n')
[err,msg] = system('wsl --status');
if err==50
    [err,msg] = system('wsl --install --no-launch');
end
if err && contains(msg,'is not recognized')
    error('%s\n%s\n%s\n%s\n',msg,...
        'Windows update 22H2 may be needed to run WSL.',...
        'If Windows update 22H2 is failing perform the install using:',...
        'https://www.microsoft.com/en-au/software-download/windows10')
elseif err
    error('%s\n',msg)
end

% Install Ubuntu
fprintf('Installing Ubuntu\n')
[err,msg] = system('wsl --install Ubuntu --no-launch');
if ~err
    [err,msg] = system('wsl --set-default Ubuntu'); %make Ubuntu the default
end
if err
    error(msg)
end

% Check BIOS
fprintf('Checking BIOS virtualization\n')
[err,msg] = system('wsl --list'); %indirect test
if err==-1 && contains(msg,'no installed distributions')
    error('%s\n%s\n%s\n%s\n','Error: CPU Virtualization may be turned off in BIOS.',...
        ' To confirm BIOS is the problem run: system("wsl --install &")',... should see: WslRegisterDistribution failed with error: 0x80370102, Please enable the Virtual Machine Platform Windows feature and ensure virtualization is enabled in the BIOS...
        ' Example HP Z8 fix: restart > F1 > Esc > Bios Setup > Security >',...
        ' > System Security > Enable "Virtualisation Technology (VTx)"')
elseif err
    error(msg)
end

% Install Astrometry.net
%In WSL install to /usr/bin/source-extractor
%Which is somwehere in C:\Program Files\WindowsApps
fprintf('%s\n','Installing astrometry.net to /usr/bin/solve-field within C:\Program Files\WindowsApps\')
[err,msg] = system('bash -c "sudo apt update"'); %may be required for next step
if ~err
    [err,msg] = system('wsl sudo apt install astrometry.net -y'); %install astrometry.net
    [~,t] = system('bash -c "dpkg -L astrometry.net | xargs file | grep executable | sed s/:.*//"'); %get all executables, wsl not working
    fprintf(' Executables: %s\n',regexprep(t,{'/usr/bin/','\n'},{'' '  '})) %list executables
    %Executables: an-fitstopnm an-pnmtofits astrometry-engine
    % build-astrometry-index downsample-fits fit-wcs fits-column-merge
    % fits-flip-endian fits-guess-scale fitsgetext get-healpix get-wcs
    % hpsplit image2xy new-wcs pad-file plot-constellations plotquad plotxy
    % query-starkd solve-field subtable tabsort wcs-grab wcs-match
    % wcs-pv2sip wcs-rd2xy wcs-resample wcs-to-tan wcs-xy2rd wcsinfo
end
if err
    error(msg)
end

% Download index files
fprintf('Downloading index files: 4200-series for %g° fields from http://data.astrometry.net\n',fov)
%The 5200-Lite may be better, its based on GAIA DR2 and Tycho2
%The 5200-Heave includes G/BP/RP mags, proper motions, parallaxes
arcmin = floor(fov/10*60); %it is recommended to get 10% of the FOV, so 5 arcmin for 0.86° fov
if ~err && arcmin <= 19
    [err,msg] = system('wsl sudo apt install astrometry.net astrometry-data-2mass-08-19 -y'); %156MB
end
if ~err && arcmin <= 7
    [err,msg] = system('wsl sudo apt install astrometry.net astrometry-data-2mass-07 -y'); %161MB
end
if ~err && arcmin <= 6
    [err,msg] = system('wsl sudo apt install astrometry.net astrometry-data-2mass-06 -y'); %328MB
end
if ~err && arcmin <= 5
    [err,msg] = system('wsl sudo apt install astrometry.net astrometry-data-2mass-05 -y'); %659MB
end
if err
    error(msg)
end

% Install Source-Extractor
fprintf('Installing Source-Extractor\n')
[err,msg] = system('wsl sudo apt install sextractor');
if err
    error(msg)
end

fprintf('Done\n')