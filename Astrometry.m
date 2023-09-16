%This class can plate solve astrophotography images by running localy
%installed <a href=http://astrometry.net>astrometry.net</a> software.
%
%Setup:
%-Download this code from: https://github.com/serg3y/matlab-astrometry
%-Add main folder to MatLab path: addpath <path>/matlab-astrometry/
%-Run: astrometry().setup
%
%-Test2:
% Astrometry().solve('D:\MatLab\matlabtoolbox.git\data\DSTG_VIS_500\32711_NAVSTAR_62_USA_201\20201214_143228.fit')
%
%Basic usage:
% as = Astrometry %create a solver, but do not solve
% as.solve(file)  %solve astrophotography image
% as.plot         %plot result.
% Additional arguments may include name/value pairs
%    (see example below):
%      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
%      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
%      radius:  approximate field size     (in deg)
%      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
%      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
%      object:  name of an object on field (string, e.g. 'M 16')
%
%Other usage:
%  as = as.LOAD(dir); as.plot;
%    Read an existing Astrometry.net set of files stored in a given directory.
%    The directory may contain WCS, CORR, RDLS, JSON, and image.
%    Then plot the result. This allows to get previous data files, or obtained
%    externally, and label them. The 'as' astrometry object must have been used
%    to solve or import astrometry data.
%  [x,y] = as.sky2xy(ra,dec)
%    Convert a RA/DEC set of coordinates (in [deg] or 'hh:mm:ss'/'deg::mm:ss')
%    into pixel coordinates on the image. The 'as' astrometry object must have
%    been used to solve or import astrometry data.
%  [ra, dec] = as.xy2sky(x,y)
%  [ra, dec] = as.xy2sky(x,y,true)
%    Convert pixel coordinates on the image into a RA/DEC set of coordinates
%    (in [deg]). When given a true argument, the result is given in
%    'hh:mm:ss'/'deg:mm:ss'. The 'as' astrometry object must have been used
%    to solve or import astrometry data.
%  f = as.findobj('object name')
%    Return information about a named object (star, deep sky object) from the
%    data base. Example: astrometry.findobj('M33')
%  as.solve(file, ...)
%    Explicitly use the local astrometry.net installation.
%    See above for the additional arguments.
%
%Using results:
% - Once an image has been solved with the 'as' object, you can use the
%   astrometry results.
% - The annotation is done asynchronously, and the Matlab prompt is recovered.
%   You may use getstatus(as) to inquire for the solve-plate status
%   (running, success, failed).
%   To wait for the end of the annotation, use as.wait. as.stop aborts it.
% - as.result.RA and as.result.Dec provide the center coordinates of the
%   field (in [deg]), while as.result.RA_hms and as.result.Dec_dms provide the
%   'HH:MM:SS' and 'Deg:MM:SS' coordinates.
% - The field rotation wrt sky is stored in as.result.rotation.
% - The pixel scale is given in [arcmin/pixel] as as.result.pixel_scale.
% - The field extension is given with its bounds as as.result.RA_min, as.result.RA_max,
%   as.result.Dec_min, and as.result.Dec_min.
% - The constellation name is stored in as.result.Constellation.
%
%Improving plate-solve efficiency:
% To facilitate the plate-solve/annotation of images, you may:
% - specify the field size with additional arguments such as:
%    astrometry(..., 'scale-low', 0.5, 'scale-high',2)
% - provide an initial guess for the location, and its range, such as:
%    astrometry('examples/M13-2018-05-19.jpg', ...
%      'ra','16:33:51','dec','30:39:35','radius', 2)
% - provide the name of an object on field, such as:
%    astrometry('examples/M13-2018-05-19.jpg','object','m 13','radius',2)
% - add more star data bases (e.g. 2MASS over Tycho2).
%
%Example:
% astrometry().solve('D:\MatLab\matlab-astrometry.git\examples\M13-2018-05-19.jpg','scale-low',0.5,'scale-high',2,'scale-high',0.9,'scale-low',0.8)
%
%Methods:
%  findobj   Find a given object in catalogs.
%  getstatus Return the astrometry status (success, failed).
%  image     Show the solve-plate image with annotations.
%  load      Load astrometry files (WCS,FITS) from a directory.
%  local     Loads an image and identifies its objects using local solve-field
%  plot      Show the solve-plate image with annotations. Same as image.
%  sky2xy    Convert RA,Dec coordinates to x,y pixels on image
%  solve     Solve an image field. Plot further results with IMAGE method.
%  stop      Ends any current annotation and reset the object.
%  visible   Return/display all visible objects on image.
%  wait      Waits for completion of the annotation.
%  xy2sky    Convert pixel image coordinates to RA,Dec.
%
%Credit:
%    sky2xy and xy2sky from E. Ofek http://weizmann.ac.il/home/eofek/matlab/
%
%(c) E. Farhi, 2018. GPL2.

%TODO:
%#ok<*TNOW1,*DATST,*TRYNC>

classdef Astrometry < handle

    properties
        file      = '' %input image file to process
        autoplot  = 0  %if true, display annotated images
        args      = [] %astrometry.net arguments for repeted use
        data_fold = '' %folder with additional star catalogs
        overwrite = 1  %if output exists overwrite it or skip it
        blocking  = 0  %wait for file to finish before resuming
        result   = []  %results from the annotation
        Timer
    end

    properties (SetAccess = protected) %read only
        status   = 'init' %can be: running, failed, success
        duration = duration();
    end

    properties (Access = private)
        process_exe = []
        process_dir = []
        timer       = []
        starttime   = []
    end

    properties (Constant = true)
        %executables = find_executables  %find installed programs
        solve_field = 'wsl solve-field'
        sextracotr %improves performance
        wcs2kml
    end

    methods
        function obj = Astrometry(varargin)
            %Initialise the class.
            % A = Astrometry()          -use default arguments
            % A = Astrometry(par1=val1,...)   -assign arguments
            %Example:
            % Astrometry().solve('M33.jpg','scale-low',0.5,'scale-high',2)

            % Parse inputs
            for k = 1:2:nargin
                obj.(varargin{k}) = varargin{k+1};
            end
            if ispc
                setenv('WSL_UTF8','1')
                %Without this system('wsl') return nulls
                %which causes output text to not display
                %https://github.com/microsoft/WSL/issues/4607#issuecomment-1197258447
                %to restore use: setenv('WSL_UTF8','')
            end
        end

        function setup(~,fov)
            %Auto install WSL, Ubuntu, Astrometry.net for Windows 10/11.
            % astrometry().setup
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
            % !wsl  -e  solve-field /mnt/c/MatLab/matlabtoolbox.git/data/DSTG_VIS_500/32711_NAVSTAR_62_USA_201/20201214_143228.fit --overwrite --downsample 2
            % !wsl      solve-field /mnt/c/MatLab/matlabtoolbox.git/data/DSTG_VIS_500/32711_NAVSTAR_62_USA_201/20201214_143228.fit --overwrite --downsample 2
            %-Ref:
            % https://www.hnsky.org/linux_subsyst.htm
            % https://github.com/Jusas/astrometry-api-lite#installation
            % https://learn.microsoft.com/en-us/windows/wsl/install
            % https://packages.debian.org/source/bookworm/astrometry-data-2mass
            % http://data.astrometry.net/

            if nargin<2 || ismept(fov), fov = 0.86; end

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
            fprintf('%s\n','Installing astrometry.net to /usr/bin/solve-field within C:\Program Files\WindowsApps\')
            [err,msg] = system('bash -c "sudo apt update"'); %may be required for next step
            if ~err
                [err,msg] = system('bash -c "sudo apt install astrometry.net -y"'); %install astrometry.net
                [~,t] = system('bash -c "dpkg -L astrometry.net | xargs file | grep executable | sed s/:.*//"'); %get all executables
                fprintf(' Executables: %s\n',regexprep(t,{'/usr/bin/','\n'},{'' '  '})) %list executables
                % an-fitstopnm  an-pnmtofits  astrometry-engine
                % build-astrometry-index  downsample-fits  fit-wcs
                % fits-column-merge  fits-flip-endian  fits-guess-scale
                % fitsgetext  get-healpix  get-wcs  hpsplit  image2xy
                % new-wcs  pad-file  plot-constellations  plotquad  plotxy 
                % query-starkd  solve-field  subtable  tabsort  wcs-grab
                % wcs-match  wcs-pv2sip  wcs-rd2xy  wcs-resample
                % wcs-to-tan wcs-xy2rd  wcsinfo
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
            [err,msg] = system('bash -c "sudo apt install sextractor"');
            if err
                error(msg)
            end

            fprintf('Done\n')
        end

        function [out,file] = solve(obj,file,varargin)
            %Solve an image field.
            % as.solve(file)      image file to process
            % as.solve(file,arg)    name/value pair, eg
            %   ra: approximate RA coordinate (deg or 'hh:mm:ss')
            %   dec: approximate DEC coordinate (deg or 'deg:mm:ss')
            %   radius: approximate field size (deg)
            %   scale-high: upper limit of the field of view (deg)
            %   scale-low:  lower limit of the field of view (deg)
            %   object: name of an object on field (string, eg 'M 16')
            %Example:
            % as.solve('M33.jpg')
            % as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','radius',2)
            % as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','scale-low',0.5,'scale-high',2)

            % Defaults
            if nargin<2 || isempty(file), file = obj.file; end
               file = convertStringsToChars(file);
            out = []; %init output

            % Check status
            if obj.Timer.Running=="on"
                disp(' Current solver is RUNNING.')
                return
            end

            % Write image data to file
            if isnumeric(file) && ~isempty(file)
                imwrite(file,'temp.png');
                file = tempname;
            end

            % Batch mode
            if isfolder(file)
                file = dir(fullfile(file,'*.fit')); %list fit files in folder as struct
            elseif ischar(file) && any(file=='*')
                file = dir([file '.fit']);  %find files using wildcard as struct
            end
            if isstruct(file)
                file = fullfile({file.folder},{file.name})'; %convert to cell
            end
            if iscell(file)
                out = cell(size(file));
                for k = 1:numel(file)
                    out{k} = obj.solve(file{k}, varargin{:});
                end
                return %finished batch mode
            end

            % Ask user to select file
            if isempty(file)
                [file, pathname] = uigetfile( ...
                    {'*.jpg;*.jpeg;*.gif;*.png;*.fit;*.fits;*.fts' 'Any image (JPG,PNG,GIF,FITS)'
                    '*.fit;*.fits;*.fts' 'FITS image (FITS)'
                    '*.jpg;*.jpeg' 'JPEG image (JPG)'
                    '*.png' 'PNG image (PNG)'
                    '*.gif' 'GIF image (GIF)'
                    '*.*', 'All Files'}, ...
                    ' Select astrophotography image');
                if isequal(file,0)
                    file = '';
                    return
                end
                file = fullfile(pathname, file);
            end

            % Resolve path
            f = dir(file);
            if isempty(f)
                error('File not found.')
            end
            file = fullfile(f.folder,f.name);
            obj.file = file;

            % Prepare folder
            if isempty(obj.process_dir) %first call
                obj.process_dir = tempname; %create temp folder name
            elseif ~isempty(dir(obj.process_dir))
                rmdir(obj.process_dir, 's'); %clean existing folder
            end
            if ~isfolder(obj.process_dir)
                mkdir(obj.process_dir); %create folder
            end
            fold = obj.process_dir;

            % Convert to unix path
            if 1 %HACK
                file = regexprep(file,'C:\','mnt/c/');
                file = regexprep(file,'\','/');
            end

            % Build the command
            cmd = [obj.solve_field];
            cmd = [cmd ' ' file];
            if ~isempty(obj.sextractor)
                cmd = [cmd ' --use-sextractor']; % highly improves annotation efficiency
            end
            cmd = [cmd ' --dir '      fold];
            cmd = [cmd ' --new-fits ' fullfile(fold, 'results.fits')];
            cmd = [cmd ' --rdls '     fullfile(fold, 'results.rdls')];
            cmd = [cmd ' --corr '     fullfile(fold, 'results.corr') ' --tag-all'];
            if ~isempty(obj.wcs2kml)
                cmd = [cmd ' --kmz '      fullfile(fold, 'results.kml') ' --no-tweak'];
            end
            cmd = [cmd ' --wcs='      fullfile(fold, 'results.wcs')];

            % Handle arguments
            if nargin > 3
                obj = [];
                rmarg = [];
                for f = 1:numel(varargin)
                    if ischar(varargin{f}) && f < numel(varargin) && any(strcmp(varargin{f}, {'findobj' 'goto' 'obj' 'object'}))
                        if ischar(varargin{f+1})
                            obj = obj.findobj(varargin{f+1});
                        elseif isstruct(varargin{f+1})
                            obj = varargin{f+1};
                        end
                        rmarg = [f f+1];
                    elseif isstruct(varargin{f})
                        obj = varargin{f};
                        rmarg = f;
                    end
                end
                if ~isempty(obj) && isstruct(obj) && isfield(obj, 'RA') && isfield(obj, 'DEC')
                    varargin(rmarg) = [];
                    varargin{end+1} = 'ra';  varargin{end+1} = obj.RA;
                    varargin{end+1} = 'dec'; varargin{end+1} = obj.DEC;
                    if isfield(obj,'NAME')
                        disp([' goto "' obj.NAME '"'])
                    end
                end
            end

            % handle additional arguments in name/value pairs
            if nargin > 3 && mod(numel(varargin), 2) == 0
                for f = 1:2:numel(varargin)
                    if ischar(varargin{f})
                        if strcmp(varargin{f}, 'scale-lower')
                            varargin{f} = 'scale-low';
                        elseif strcmp(varargin{f}, 'scale-upper')
                            varargin{f} = 'scale-high';
                        end
                        cmd = [cmd ' --' varargin{f} '=' num2str(varargin{f+1})];
                    end
                end
            end

            % Execute command
            obj.Timer = startcmd(cmd);
        end

        function ret = load(obj, d, varargin)
            %Load astrometry files (WCS,FITS) from a directory
            %  LOAD(astrometry, directory);
            %  The directory may contain WCS, CORR, RDLS or JSON, and image.
            %  No solve plate is performed, only data is read.
            %
            %  LOAD(astrometry, image, ...) starts solve-plate annotation, same as SOLVE.
            if nargin < 2
                d = obj.process_dir;
            end
            if ~isfolder(d)
                % loading an image ?
                try
                    % im = imread(d);
                    solve(obj, d, varargin{:});
                    return
                end
            end

            obj.result = getresult(d, obj);
            if isempty(obj.result)
                obj.status = 'failed';
                obj.process_exe = [];
            else
                obj.status = 'success';
                % is the image available ? use one from the directory
                if ~exist(obj.file, 'file')
                    % search in the result directory
                    d = [dir(fullfile(obj.result.dir, '*.png')); dir(fullfile(obj.result.dir, '*.fits'))];
                    if ~isempty(d)
                        d = d(1);
                        obj.file = fullfile(obj.result.dir, d.name);
                    end
                end
            end
            ret = obj.result;
        end

        function ret = getstatus(obj)
            %Return astrometry status, ie 'success' 'failed'
            ret = obj.status;
        end

        function fig = plot(obj, mag)
            %Show the solve-plate image with annotations
            %  as.plot
            %  as.plot(mag)   limits the objects up to given magnitude
            %Example:
            %  astrometry(file).plot

            if nargin < 2, mag = inf; end

            fig = [];
            if ~ischar(obj.file) || isempty(obj.file) || isempty(dir(obj.file))
                return
            end
            try
                im  = imread(obj.file);
            catch ex
                getReport(ex)
                disp([mfilename  ': ERROR: can not read image ' obj.file])
                return
            end
            fig = figure('Name', [' ' obj.file]);
            image(im)
            clear im

            ret = obj.result;

            % Set title
            [~,f,e] = fileparts(obj.file);
            if isfield(ret, 'Constellation')
                title([f e ' in ' ret.Constellation]);
            else
                title([f e]);
            end

            % Overlay results
            if ~isempty(obj.result) && strcmp(obj.status, 'success') && isfield(obj.result, 'RA_hms')

                hold on
                % Central coordinates
                sz = obj.result.size/2;
                h  = plot(sz(1), sz(2), 'r+'); set(h, 'MarkerSize', 16);
                hcmenu = uicontextmenu;
                uimenu(hcmenu, 'Label', '<html><b>Field center</b></html>');
                uimenu(hcmenu, 'Label', ['RA=  ' ret.RA_hms]);
                uimenu(hcmenu, 'Label', ['DEC= ' ret.Dec_dms]);
                uimenu(hcmenu, 'Label', ['Rotation= ' num2str(ret.rotation) ' [deg]']);
                set(h, 'UIContextMenu', hcmenu);

                % Get list of visible objects
                v = visible(obj, mag);

                % Make sure we do not plot too many
                if numel(v) > 1500
                    v = v(1:1500);
                end

                for index = 1:numel(v)
                    % Find all objects from data base within bounds
                    this = v(index);

                    % Stars in green, DSO in cyan
                    if strcmp(this.catalog,'stars')
                        c = 'g';
                    else
                        c = 'c';
                    end
                    x = this.X;
                    y = this.Y;

                    % Plot symbol
                    if isfinite(this.SIZE) && this.SIZE > 5
                        h = plot(x, y, [c 'o']);
                        set(h, 'MarkerSize', ceil(this.SIZE));
                    else
                        h = plot(x, y, [c 's']);
                        this.SIZE=12;
                    end

                    % Context menu
                    hcmenu = uicontextmenu;
                    uimenu(hcmenu, 'Label', ['RA=  ' this.RA ' [' num2str(this.RA_deg) ' deg]']);
                    uimenu(hcmenu, 'Label', ['DEC= ' this.DEC ' [' num2str(this.DEC_deg) ' deg]']);
                    uimenu(hcmenu, 'Label', ['<html><b>' this.NAME '</html></b>'], 'Separator', 'on');
                    uimenu(hcmenu, 'Label', ['TYPE: ' this.TYPE]);
                    if isfinite(this.MAG) && this.MAG > 0
                        uimenu(hcmenu, 'Label', ['MAGNITUDE= ' num2str(this.MAG)]);
                    end
                    if isfinite(this.DIST) && this.DIST > 0
                        uimenu(hcmenu, 'Label', ['DIST= ' sprintf('%.3g', this.DIST*3.262) ' [ly]']);
                    end
                    set(h, 'UIContextMenu', hcmenu);
                    if numel(v) < 1500
                        text(x+this.SIZE,y-this.SIZE,this.NAME,'Color', c);
                    end
                end
            end
        end

        function plot2(obj)
            % Display Astrometry object (short)

            if ~isempty(inputname(1))
                iname = inputname(1);
            else
                iname = 'ans';
            end
            if isdeployed || ~usejava('jvm') || ~usejava('desktop')
                id = class(obj);
            else
                id = ['<a href="matlab:doc ' class(obj) '">' class(obj) '</a> ' ...
                    '(<a href="matlab:methods ' class(obj) '">methods</a>,' ...
                    '<a href="matlab:image(' iname ');">plot</a>,' ...
                    '<a href="matlab:disp(' iname ');">more...</a>)'];
            end
            if obj.Timer.Running=="on"
                fprintf('%s = %s for "%s" BUSY\n',iname, id, obj.file)
            else
                fprintf('%s = %s for "%s"\n',iname, id, obj.file)
            end
        end

        function [x,y] = sky2xy(obj, ra, dec)
            % Convert RA,Dec (deg) to x,y image pixel coordinates.
            x = [];
            y = [];
            if isempty(obj.result)
                return
            end
            if ~isscalar(ra)
                ra = getra(ra);
            end
            if ~isscalar(dec)
                dec = getdec(dec);
            end
            [x,y] = sky2xy_tan(obj.result.wcs.meta, ra*pi/180, dec*pi/180); % MAAT Ofek (private)
        end

        function [ra,dec] = xy2sky(obj, x, y, str)
            % Convert x,y image pixel to RA,Dec (deg) coordinates.
            ra = [];
            dec = [];
            if isempty(obj.result)
                return
            end
            if nargin > 3 && str == "string"
                str = true;
            else
                str = false;
            end
            [ra, dec] = xy2sky_tan(obj.result.wcs.meta, x,y); % MAAT Ofek (private)
            ra  = rad2deg(ra);
            dec = rad2deg(dec);
            if str
                ra = getra (ra/15, true);
                dec= getdec(dec,   true);
            end
        end

        function out = findobj(obj,name)
            %Search star catalogs for an object by name
            % out = as.findobj(name)
            out = []; %init
            for f = obj.catalogs'
                data = obj.catalogs(f{1});
                if isfield(data,'MAG')
                    ind = find(contains(regexprep(data.NAME,' ',''),name,'IgnoreCase',true)); %search for name
                    if ~isempty(ind)
                        out.catalog = f{1};
                        out.index   = ind';
                        out.NAME    = data.NAME(ind)';
                        out.RA      = data.RA(ind)';
                        out.DEC     = data.DEC(ind)';
                        out.MAG     = data.MAG(ind)';
                        out.TYPE    = data.TYPE(ind)';
                        out.DIST    = data.DIST(ind)';
                        break
                    end
                end
            end
            % Print results
            if ~nargout
                if ~isempty(out)
                    fprintf('Found %g matches for "%s":\n',numel(ind),name)
                    disp(out)
                else
                    fprintf('Filed to find "%s"\n',name)
                end
                clear out
            end
        end

        function v = visible(obj, mag)
            %List visible objects in image
            %  visible(as)
            %  visible(as,mag) limits the objects up to given magnitude

            v = [];
            if nargin < 2, mag = inf; end

            if ~isempty(obj.result) && isfield(obj.result, 'RA_hms')

                ret = obj.result;

                for catalogs_names = {'stars' 'deep_sky_objects'}
                    % find all objects from data base within bounds
                    catalog = obj.catalogs.(catalogs_names{1});

                    found = find(ret.RA_min <= catalog.RA & catalog.RA <= ret.RA_max & ret.Dec_min<= catalog.DEC & catalog.DEC <= ret.Dec_max & catalog.MAG <= mag);

                    if ~isempty(found)
                        disp([mfilename  ': ' num2str(numel(found)) ' referenced ' catalogs_names{1} ' in field.'])
                    end
                    for index = 1:numel(found)
                        obj   = found(index);
                        ra    = catalog.RA(obj);
                        dec   = catalog.DEC(obj);
                        [x,y] = obj.sky2xy(ra, dec);

                        % ignore when not on image
                        if x < 1 || x > ret.size(1) || y < 1 || y > ret.size(2)
                            continue
                        end

                        this.RA   = getra(ra/15, true);
                        this.DEC  = getdec(dec,  true);
                        this.NAME = catalog.NAME{obj};
                        this.TYPE = catalog.TYPE{obj};
                        this.MAG  = catalog.MAG(obj);
                        this.SIZE = catalog.SIZE(obj); % arcmin
                        this.DIST = catalog.DIST(obj);
                        this.X    = x;
                        this.Y    = y;
                        this.catalog = catalogs_names{1};
                        this.RA_deg  = ra;
                        this.DEC_deg = dec;

                        if isempty(v)
                            v = this;
                        else
                            v(end+1) = this;
                        end
                    end
                end

                % sort by increasing magnitude (decreasing brightness)
                [~,index]=sort([v.MAG]);
                v = v(index);

                % display the list
                if nargout == 0
                    disp(obj.file)
                    disp 'TYPE            MAG  RA              DEC                 DIST  NAME'
                    for index = 1:numel(v)
                        this = v(index);
                        fprintf('%-12s  %5.1f  %-14s  %-14s  %8.2g  %s\n', this.TYPE, this.MAG, this.RA, this.DEC, this.DIST*3.262, this.NAME) % display the list
                    end
                end
            end
        end

        function disp(obj)
            % DISP Display Astrometry object (details)

            if ~isempty(inputname(1))
                iname = inputname(1);
            else
                iname = 'ans';
            end
            if isdeployed || ~usejava('jvm') || ~usejava('desktop'), id=class(obj);
            else
                id = ['<a href="matlab:doc ' class(obj) '">' class(obj) '</a> ' ...
                    '(<a href="matlab:methods ' class(obj) '">methods</a>,' ...
                    '<a href="matlab:image(' iname ');">plot</a>,' ...
                    '<a href="matlab:visible(' iname ');">visible...</a>)'];
            end
            fprintf('%s = %s for "%s":\n', iname, id, obj.file)
            if ~isempty(obj.result) && strcmp(obj.status, 'success') && isfield(obj.result, 'RA_hms')
                % get list of visible objects, extract brightest star and DSO
                v = visible(obj);
                mag = [inf inf];  % star, dso
                name= {''  ''};

                for index = 1:numel(v)
                    % find all objects from data base within bounds
                    this = v(index);
                    % stars in green, DSO in cyan
                    if strcmp(this.catalog, 'stars')
                        obj_index = 1;
                    else
                        obj_index = 2;
                    end
                    if this.MAG < mag(obj_index)
                        mag(obj_index)  = this.MAG;
                        name{obj_index} = [this.NAME ' ' this.TYPE];
                    end
                end
                if isfield(obj.result, 'Constellation')
                    disp(['  Constellation: ' obj.result.Constellation])
                end
                disp(['  RA:            ' obj.result.RA_hms  ' [h:min:s]; '   num2str(obj.result.RA)  ' [deg]'])
                disp(['  DEC:           ' obj.result.Dec_dms ' [deg:min:s]; ' num2str(obj.result.Dec) ' [deg]'])
                for index = 1:2
                    if ~isempty(name{index})
                        disp(['    magnitude ' num2str(mag(index)) ' ' name{index}])
                    end
                end
                disp(['  Rotation:      ' num2str(obj.result.rotation)    ' [deg] (to get sky view)']);
                disp(['  Pixel scale:   ' num2str(obj.result.pixel_scale) ' [arcsec/pixel]']);
                if isdeployed || ~usejava('jvm') || ~usejava('desktop')
                    disp(['  Results are in ' obj.process_dir]);
                else
                    disp(['  Results are in <a href="' obj.process_dir '">' obj.process_dir '</a>']);
                end
                builtin('disp',obj)
                disp([iname '.result:'])
                disp(obj.result);
            elseif ~isempty(obj.process_dir) && isfolder(obj.process_dir)
                if isdeployed || ~usejava('jvm') || ~usejava('desktop')
                    disp(['  ' upper(obj.status) ' in ' obj.process_dir]);
                else
                    disp(['  ' upper(obj.status) ' in <a href="' obj.process_dir '">' obj.process_dir '</a>']);
                end
            else
                disp(['  ' upper(obj.status) ': use as.annotate(''file'').']);
            end

        end

        function wait(obj)
            % Waits for solve
            while obj.Timer.Running=="on"
                pause(1)
            end
        end

        function stop(obj)
            % Stop current solve
            if ~isempty(obj.Timer)
                obj.Timer.stop
                obj.result = obj.Timer.UserData;
                obj.Timer.delete
            end
        end

        function [out,paths] = catalogs(obj,name)
            %Load star catalog data
            % [list,paths] = as.catalogs  -list catalog names & file paths
            % data = as.catalogs(name)
            %Example:
            % list = astrometry().catalogs
            % data = astrometry().catalogs('constellations')
            % data = astrometry().catalogs(astrometry().catalogs{2})
            persistent old_data
            if isempty(obj.data_fold)
                obj.data_fold = fullfile(fileparts(mfilename('fullpath')),'catalogs');
            end
            if nargin<2 || isempty(name)
                paths = dir(fullfile(obj.data_fold,'*.mat')); %find mat files
                [~,out] = fileparts({paths.name}'); %file names
                paths = fullfile({paths.folder},{paths.name})'; %file paths
            elseif ~isempty(old_data) && isfield(old_data,name)
                out = old_data.(name);
            else
                out = load(fullfile(obj.data_fold,name)).(name); %load catalog
                if ~isfield(old_data,name)
                    old_data.(name) = out; %cache for next time
                end
            end
        end

    end

end


function out = getresult(d, obj)
% getresult: extract WCS and star matching information from the output files.
%
% input:
%   d: directory where astrometry.net results are stored.

out = [];
for file = {'results.wcs' 'wcs.fits'}
    if exist(fullfile(d, file{1}), 'file')
        out.wcs  = read_fits(fullfile(d, file{1}));

        % get image center and print it
        if isfield(out.wcs,'meta') && isfield(out.wcs.meta,'CRVAL1')
            wcs = out.wcs.meta;
            out.wcs.meta.CD = [wcs.CD1_1 wcs.CD1_2 ; wcs.CD2_1 wcs.CD2_2];

            % get central coordinates
            out.size = [wcs.IMAGEW wcs.IMAGEH];
            sz  = out.size/2;

            [out.RA, out.Dec] = xy2sky_tan(out.wcs.meta, sz(1), sz(2)); % MAAT Ofek (private)
            out.RA       = rad2deg(out.RA);
            out.Dec      = rad2deg(out.Dec);
            out.RA_hms   = getra(out.RA/15,true);
            out.Dec_dms  = getdec(out.Dec, true);

            out.pixel_scale = sqrt(abs(wcs.CD1_1 * wcs.CD2_2  - wcs.CD1_2 * wcs.CD2_1))*3600; %pixel scale (arcsec/pixel)
            out.rotation = atan2d(wcs.CD2_1, wcs.CD1_1); %rotation angle

            % compute RA,Dec image bounds
            RA = []; Dec = [];
            [RA(end+1), Dec(end+1)] = xy2sky_tan(out.wcs.meta, wcs.IMAGEW, wcs.IMAGEH);
            [RA(end+1), Dec(end+1)] = xy2sky_tan(out.wcs.meta, 1         , wcs.IMAGEH);
            [RA(end+1), Dec(end+1)] = xy2sky_tan(out.wcs.meta, wcs.IMAGEW, 1);
            [RA(end+1), Dec(end+1)] = xy2sky_tan(out.wcs.meta, 1         , 1);
            RA  = rad2deg(RA);
            Dec = rad2deg(Dec);
            out.RA_min  = min(RA);
            out.RA_max  = max(RA);
            out.Dec_min = min(Dec);
            out.Dec_max = max(Dec);

            %Find nearest constellation
            cat = obj.catalogs('constellations');
            [~,ind] = min((out.RA - cat.RA).^2 + (out.Dec - cat.DEC).^2);
            out.Constellation = cat.Name{ind};
        end
    end
end
for file = {'results.rdls' 'rdls.fits'}
    if exist(fullfile(d, file{1}), 'file')
        out.rdls = read_fits(fullfile(d, file{1}));
    end
end
for file = {'results.corr' 'corr.fits'}
    if exist(fullfile(d, file{1}), 'file')
        out.corr = read_fits(fullfile(d, file{1}));
    end
end
if exist(fullfile(d, 'results.json'), 'file')
    out.json = loadjson(fullfile(d, 'results.json'));
end
if ~isempty(out)
    out.dir     = d;
end
end


function executables = find_executables
% locate executables, return a structure

persistent executables_cache % stored here so that they are not searched for further calls

if ~isempty(executables_cache)
    executables = executables_cache;
else
    if     ismac,  cmd_prefix = 'DYLD_LIBRARY_PATH= ;';
    elseif isunix, cmd_prefix = 'LD_LIBRARY_PATH= ; ';
    else,          cmd_prefix = '';
    end
    if ispc, ext = '.exe';
    else,    ext = '';      %linux
    end
    root_fold = fileparts(mfilename('fullpath'));

    % Searching for installed programs
    fprintf('Searching for available programs:\n') %progress
    for exe = ["solve-field" "sextractor" "wcs2kml"] %try these programs
        for cmd = [fullfile(root_fold,exe+ext) fullfile(root_fold,exe) exe+ext exe] %try these calls
            if isfile(cmd)
                err = 0;
            else
                [err,~] = system(cmd_prefix + cmd + " --version"); %run from Matlab
            end
            name = regexprep(exe, {'-' '\.'}, '_');
            if ~err %executable is found, possible errors 1,127,9009
                executables.(name) = cmd;
                break
            else
                executables.(name) = [];
            end
        end
        fprintf(' %12s: %s\n', exe, executables.(name)) %progress
    end
    executables_cache = executables; %cache results
end
end

function [ra_h, ra_min, ra_s] = getra(ra, asstring)
% getra: convert any input RA (in hours) into h and min

if nargin < 2 || isempty(asstring), asstring = false; end
if ischar(ra)
    ra = repradec(ra);
end
if isstruct(ra) && isfield(ra, 'RA')
    ra = ra.RA;
end
if isstruct(ra)
    if isfield(ra, 'h'),   ra_h   = ra.h; end
    if isfield(ra, 'min'), ra_min = ra.min; end
    if isfield(ra, 's'),   ra_s   = ra.s; end
end
if isnumeric(ra)
    if isscalar(ra)
        ra_h   = fix(ra);
        ra_min = abs(ra - ra_h)*60;
        ra_s   = abs(ra_min - fix(ra_min))*60;
        ra_min = fix(ra_min);
    elseif numel(ra) == 2
        ra_h   = ra(1);
        ra_min = abs(ra(2));
        ra_s   = abs(ra_min - fix(ra_min))*60;
        ra_min = fix(ra_min);
    elseif numel(ra) == 3
        ra_h   = ra(1);
        ra_min = ra(2);
        ra_s   = ra(3);
    end
else
    disp([mfilename  ': invalid RA.'])
    disp(ra)
end
if nargout == 1
    if asstring
        ra_h = [num2str(ra_h) ':' num2str(ra_min) ':' num2str(ra_s)];
    else
        ra_h = ra_h+ra_min/60 + ra_s/3600;
    end
elseif nargout == 2
    ra_min = ra_min + ra_s/60;
end
end

function [dec_deg, dec_min, dec_s] = getdec(dec, asstring)
% Convert any input DEC into deg and min

if nargin < 2 || isempty(asstring), asstring = false; end
if ischar(dec)
    dec = repradec(dec);
end
if isstruct(dec) && isfield(dec, 'DEC')
    dec = dec.DEC;
end
if isstruct(dec)
    if isfield(dec, 'deg'), dec_deg = dec.deg; end
    if isfield(dec, 'min'), dec_min = dec.min; end
    if isfield(dec, 'min'), dec_s   = dec.s; end
end
if isnumeric(dec)
    if isscalar(dec)
        dec_deg = floor(dec);
        dec_min = abs(dec - dec_deg)*60;
        dec_s   = abs(dec_min - fix(dec_min))*60;
        dec_min = fix(dec_min);
    elseif numel(dec) == 2
        dec_deg = dec(1);
        dec_min = abs(dec(2));
        dec_s   = abs(dec_min - fix(dec_min))*60;
        dec_min = fix(dec_min);
    elseif numel(dec) == 3
        dec_deg = dec(1);
        dec_min = dec(2);
        dec_s   = dec(3);
    end
else
    fprintf(2,' invalid DEC: %g\n',dec)
end
if nargout == 1
    if asstring
        dec_deg = [num2str(dec_deg) ':' num2str(dec_min) ':' num2str(dec_s)];
    else
        dec_deg = dec_deg + dec_min/60 + dec_s/3600;
    end
elseif nargout == 2
    dec_min = dec_min + dec_s/60;
end
end

function str = repradec(str)
% repradec: replace string stuff and get it into num
str = lower(str);
for rep = {'h' 'm' 's' ':' '°' 'deg' 'd' '''' '"'}
    str = strrep(str, rep{1}, ' ');
end
str = str2num(str);
end

function ret = open_system_browser(url)
% opens URL with system browser. Returns non zero in case of error.
if strncmp(url, 'file://', length('file://'))
    url = url(8:end);
end
if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; ';
else,          precmd = '';
end
if ispc
    ret = system([precmd 'start "' url '"']);
elseif ismac
    ret = system([precmd 'open "' url '"']);
else
    [ret,~] = system([precmd 'xdg-open "' url '"']);
end
end

function varargout = startcmd(cmd,vrb)
%Run a system command as a separate process, which can be monitored.
% startcmd(cmd)       -start command
% startcmd(cmd,vrb)   -verbose: 0=none, 1=show progress, 2=show results
% T = startcmd(__)    -return timer for monitoring the process
%
%Remarks:
%-If output T is returned then it is up to the user to delete the timer
% when it is no longer required using T.delete
%-Use T.UserData to access process output when the process ends.
%-To stop the process prematurely use T.stop
%
%Example:
% startcmd('ping localhost',2)  %start process with verbose=2
%
%Example:
% T = startcmd('ping localhost');  %start process and return the timer
% while T.Running=="on"
%     pause(0.1)  %wait for process to finish
% end
% T.UserData  %process output
% T.delete  %delete timer

if nargin<2 || isempty(vrb), vrb = 0; end
delTimer = ~nargout; %delete timer when the process ends

Proc = java.lang.Runtime.getRuntime.exec(cmd); %start process
Timer = timer(TimerFcn=@Update,StopFcn=@Stop,ExecutionMode='fixedSpacing',Period=0.1); %timer to monitor the process
Timer.start; %start timer

if nargout
    varargout = {Timer}; %return timer, if requested
end

    function Update(~,~)
        if ~get(Proc,'Alive')
            Timer.stop %stop the timer, which also calls the Stop function
        end
    end

    function Stop(~,~)
        Proc.destroy %close java runtime
        t = java.util.Scanner(java.io.InputStreamReader(Proc.getInputStream)).useDelimiter('\A').next; %read process output
        Timer.UserData = regexprep(strtrim(char(t)),'\r',''); %save output to Timer.UserData
        if vrb >= 1, fprintf('Finished: %s\n',cmd), end %show progress
        if vrb >= 2, fprintf('%s\n',Timer.UserData), end %show output
        if delTimer, Timer.delete, end %delete timer
    end
end

% Solve the given astrophotography image.
% Then plot the result. Additional arguments may include name/value pairs
% (see example below):
%
%        - ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
%        - dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
%        - radius:  approximate field size     (in deg)
%        - scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
%        - scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
%
% These two syntaxes will try first any local astrometry.net installation, and
% if failed, the http://nova.astrometry.net/ service.
%
%    **as.plot**
%
% Plot the astrometry solution. The identified objects are indicated in green for
% stars, cyan circle for extended deep sky objects, and cyan squares for localized
% deep sky objects. Each indicated object has a contextual menu which gives more
% information (right click). The central coordinate of the field is shown in red.
%
% Going further <a id=going-further></a>
% =============
%
%    **as = as.load(dir); as.plot;**
%
% Read an existing Astrometry.net set of files stored in a given directory.
% The directory may contain WCS, CORR, RDLS, JSON, and image.
% Then plot the result. This allows to get previous data files, or obtained
% externally, and label them. The 'as' astrometry object must have been used
% to solve or import astrometry data.
%
%    **[x,y] = as.sky2xy(ra, dec)**
%
% Convert a RA/DEC set of coordinates (in [deg] or 'hh:mm:ss'/'deg::mm:ss')
% into pixel coordinates on the image. The 'as' astrometry object must have
% been used to solve or import astrometry data.
%
%    **[ra, dec] = as.xy2sky(x,y)**
%
%    **[ra, dec] = as.xy2sky(x,y, 'string')**
%
% Convert pixel coordinates on the image into a RA/DEC set of coordinates
% (in [deg]). When given a 'string' argument, the result is given in
% 'hh:mm:ss'/'deg:mm:ss'. The 'as' astrometry object must have been used
% to solve or import astrometry data.
%
%    **f = as.FINDOBJ('object name')**
%
% Return information about a named object (star, deep sky object) from the
% data base. Example: astrometry().findobj('M33')
%
%    **as.solve(file, ...);**
%
%
%
% Using results <a id=using-results></a>
% =============
% Once an image has been solved with the 'as' object, you can use the astrometry results.
%
% - **as.result.RA** and **as.result.Dec** provide the center coordinates of the
%   field (in [deg]), while **as.result.RA_hms** and **as.result.Dec_dms** provide the
%   'HH:MM:SS' and 'Deg:MM:SS' coordinates.
% - The field rotation wrt sky is stored in **as.result.rotation**.
% - The pixel scale is given in [arcmin/pixel] as **as.result.pixel_scale**.
% - The field extension is given with its bounds as **as.result.RA_min**, **as.result.RA_max**,
%   **as.result.Dec_min**, and **as.result.Dec_min**.
% - The constellation name is stored in **as.result.Constellation**.
%
% Improving the plate-solve efficiency <a id=improving-efficiency></a>
% ====================================
%
% To facilitate the plate-solve/annotation of images, you may:
%
% - specify the field size with additional arguments such as:
%   `astrometry(..., 'scale-low', 0.5, 'scale-high',2)`
%   This is what **works best**, by far.
%
% - provide an initial guess for the location, and its range, such as:
%   `astrometry('examples/M33-2018-05-19.jpg','ra','01:33:51','dec','30:39:35','radius', 2)`
%
% - provide the name of an object on field, such as:
%  `astrometry('examples/M13-2018-05-19.jpg','object','m 13','radius',2)`
%
% - add more star data bases (e.g. 2MASS over Tycho2).
%
% Examples <a id=examples></a>
% ========
%
% ```matlab
%   as=astrometry('examples/M33-2018-08-15.jpg','scale-low', 0.5, 'scale-high',2);
%   as.image; % once done
% ```
%
%   You will then get, in about 30 sec, the image:
%   ![Image of Astrometry](https://github.com/farhi/matlab-astrometry/blob/master/examples/M33-solved.png)
%
% Methods <a id=methods></a>
% - findobj   find a given object in catalogs.
% - getstatus return the astrometry status (success, failed)
% - image     show the solve-plate image with annotations
% - load      load astrometry files (WCS,FITS) from a directory
% - local     loads an image and identifies its objects using local solve-field
% - plot      show the solve-plate image with annotations. Same as image.
% - sky2xy    convert RA,Dec coordinates to x,y pixels on image
% - solve     solve an image field. Plot further results with IMAGE method.
% - stop      ends any current annotation and reset the object.
% - visible   return/display all visible objects on image
% - wait      waits for completion of the annotation
% - xy2sky    convert pixel image coordinates to RA,Dec
%
% Installation <a id=installation></a>
%    **Local installation (recommended)**
%
% On Linux systems, install the astrometry.net package, as well as the 'tycho2' data base. On Debian-class systems, this is achieved with:
%
% ```bash
%   sudo apt install astrometry.net astrometry-data-tycho2 sextractor
% ```
%
% On other systems, you will most probably need to compile it.
% See: http://astrometry.net/doc/build.html
% RedHat/Arch/MacOSX have specific installation instructions.
%
% If you have images spanning on very tiny areas (e.g. much smaller than a
% degree), you will most probably need to install the '2MASS' data base.
%
%    **Matlab files**
%
% First navigate to the matlab-astrometry directory or type:
%
% ```matlab
%   addpath /path/to/matlab-astrometry
% ```
%
% Credits <a id=credits></a>
% =======
%
% **sky2xy** and **xy2sky** from E. Ofek http://weizmann.ac.il/home/eofek/matlab/
%
% (c) E. Farhi, 2019. GPL2.

%NOTESL
% fprintf('Installing WSL, administrator privilages required...\n')
% [a,b]=system('wsl --install --distribution Ubuntu --no-launch')

%prior wsl --install ..
% setenv('WSL_UTF8','1');
% [a,b]=system('wsl') %1 13
% [a,b]=system('wsl --help') %1    13
% [a,b]=system('wsl --status') %50   []

%just after wsl --install ..
% [a,b]=system('wsl') %-1
% b =
%     'This application requires the Windows Subsystem for Linux Optional Component.
%      The system may need to be restarted so the changes can take effect.
%      Error code: Wsl/WSL_E_WSL_OPTIONAL_COMPONENT_REQUIRED
%      '
% [a,b]=system('wsl --help') %-1 lots of text, n=5667
% [a,b]=system('wsl --status') %-1
% b =
%     'This application requires the Windows Subsystem for Linux Optional Component.
%      The system may need to be restarted so the changes can take effect.
%      Error code: Wsl/WSL_E_WSL_OPTIONAL_COMPONENT_REQUIRED
%      '

%after reboot
% [a,b]=system('wsl --install --no-distribution') % -1, lots of text n=5667