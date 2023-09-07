%This class is a wrapper for locally installed http://astrometry.net
%software or the web based API, used for solve astrophotography images.
%
%Setup on Linux:
% -On Debian/Ubuntu Linux run:
%  sudo apt install astrometry.net astrometry-data-tycho2 sextractor
% -To solve fields smaller than a degree also install '2MASS' database:
%  ??? sudo apt install ???
% -Other Linux destributions (RedHat/Arch/MacOSX) require compilation:
%  See http://astrometry.net/doc/build.html
%
%Setup WSL on windows:
% See https://www.hnsky.org/linux_subsyst.htm
% See https://learn.microsoft.com/en-us/windows/wsl/install
% 1) Install WSL Ubuntu using windows command line, then REBOOT COMPUTER:
%     wsl --install         %required 1st time for other commands to work!
%     wsl --status          %show if wsl is installed, blank otherwise
%     wsl --list            %list already installed versions
%     wsl --list --online   %list versions that can be installed
%     wsl --install Ubuntu -n  %install Ubuntu but do not launch it
%     wsl                   %run the default version or print help
%     wsl --install -d Ubuntu
% 2) Test run a Linux command from windows cmd:
%     wsl -e ls   
%     wsl -- ls
%     wsl ls
%
%Setup Cygwin on Windows
% WRITE ME
%
%Web setup
%   1) Install Python
%   2) Create a 'NOVA astrometry API' key
%   3) Enter the KEY when prompted, or set it with:
%    as = astrometry;
%    as.api_key = 'mykey';
%    as.web(file, ...)
%
%Setup:
% 1) Download or clone matlab-astrometry from github:
%    https://github.com/serg3y/matlab-astrometry
% 2) Add rootfolder to MatLab path:
%    addpath <path>/matlab-astrometry
%
%Basic usage:
%  as = astrometry;
%    Create a solver, but does not solve.
%    Use solve(as, file) or local(as, file) or web(as, file) afterwards.
%  as = astrometry(file, ...); image(as);
%    Solve the given astrophotography image with local or web method.
%    Then plot the result. Additional arguments may include name/value pairs
%    (see example below):
%      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
%      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
%      radius:  approximate field size     (in deg)
%      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
%      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
%      object:  name of an object on field (string, e.g. 'M 16')
% - These two syntaxes will try first any local astrometry.net installation,
%   and if failed, the http://nova.astrometry.net/ service.
% - When the annotation has ended, an 'annotationEnd' event is triggered. You
%   may monitor this with e.g.:
%    as = astrometry('examples/M13-2018-05-19.jpg');
%    addlistener(as, 'annotationEnd', @(src,evt)disp('annotation just end'))
%
%Other usage:
%  as = LOAD(as, dir); IMAGE(as);
%    Read an existing Astrometry.net set of files stored in a given directory.
%    The directory may contain WCS, CORR, RDLS, JSON, and image.
%    Then plot the result. This allows to get previous data files, or obtained
%    externally, and label them. The 'as' astrometry object must have been used
%    to solve or import astrometry data.
%  [x,y] = SKY2XY(as, ra, dec)
%    Convert a RA/DEC set of coordinates (in [deg] or 'hh:mm:ss'/'deg::mm:ss')
%    into pixel coordinates on the image. The 'as' astrometry object must have
%    been used to solve or import astrometry data.
%  [ra, dec] = XY2SKY(as, x,y)
%  [ra, dec] = XY2SKY(as, x,y, true)
%    Convert pixel coordinates on the image into a RA/DEC set of coordinates
%    (in [deg]). When given a true argument, the result is given in
%    'hh:mm:ss'/'deg:mm:ss'. The 'as' astrometry object must have been used
%    to solve or import astrometry data.
%  f = FINDOBJ(as,'object name')
%    Return information about a named object (star, deep sky object) from the
%    data base. Example: astrometry.findobj('M33')
%  LOCAL(as, file, ...)
%    Explicitly use the local 'solve-field' astrometry.net installation.
%    See above for the additional arguments.
%  WEB(as, file, ...)
%    Explicitly use the http://nova.astrometry.net/ web service.
%    See above for the additional arguments.
%  WEB(as)
%    For a solved image, the corresponding sky view is displayed on
%    http://www.sky-map.org . The 'as' astrometry object must have been used
%    to solve or import astrometry data.
%
%Using results:
% - Once an image has been solved with the 'as' object, you can use the
%   astrometry results.
% - The annotation is done asynchronously, and the Matlab prompt is recovered.
%   You may use getstatus(as) to inquire for the solve-plate status
%   (running, success, failed).
%   To wait for the end of the annotation, use waitfor(as). stop(as) aborts it.
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
% as=astrometry('examples/M13-2018-05-19.jpg','scale-low', 0.5, 'scale-high',2);
% image(as);
%
%Methods
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
% - waitfor   waits for completion of the annotation
% - web       loads an image and identifies its objects using web service
% - xy2sky    convert pixel image coordinates to RA,Dec
%
%Credit:
%    sky2xy and xy2sky from E. Ofek http://weizmann.ac.il/home/eofek/matlab/
%
%(c) E. Farhi, 2018. GPL2.

%TODO:
%#ok<*TNOW1,*DATST,*TRYNC>

classdef astrometry < handle

    properties
        api_key    = ''     %api-key for nova.astrometry.net
        % example: 'kvfubnepntofzpcl' 'ghqpqhztzychczjh'
        % from: https://git.kpi.fei.tuke.sk/TP/ExplorationOfInterstellarObjects/blob/master/src/sk/tuke/fei/kpi/tp/eoio/AstrometryAPI.java
        result     = []     %results from the annotation or empty when failed
        filename   = ''     %the image to annotate
        status     = 'init' %can be: running, failed, success
        catalogs   = []     %catalogs of common objects
        vargin     = []     %arguments stored at instantiation for reuse
        autoplot   = false  %when true, display annotated image on success
        duration   = 0
    end

    properties (Access = private)
        process_java = []
        process_dir  = []
        timer        = []
        starttime    = []
    end

    properties (Constant = true)
        executables = find_executables  %find installed programs
    end

    events
        annotationStart
        annotationEnd
        idle
        busy
    end

    methods
        function obj = astrometry(filename, varargin)
            % ASTROMETRY Loads an image and identifies its objects using astrometry.net
            %
            % as = astrometry;
            %   Create a solver, but does not solve.
            %   Use local(as, file) or web(as, file) afterwards
            % as = astrometry(file, ...);
            %   Solve the given astrophotography image with local or web method.
            %
            % input(optional):
            %    filename: an image to annotate
            %    any name/value pair as:
            %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
            %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
            %      radius:  approximate field size     (in deg)
            %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
            %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
            %      object:  name of an object on field (string, e.g. 'M 16')
            %
            % Example:
            %   as=astrometry('M33.jpg','scale-low', 0.5, 'scale-high',2);

            % handle input arguments
            removeme = [];
            for index = 1:2:numel(varargin)
                if ischar(varargin{index})
                    switch varargin{index}
                        case 'catalogs'
                            obj.catalogs = varargin{index+1};
                            removeme = [removeme index index+1];
                        case 'autoplot'
                            obj.autoplot = true;
                            removeme = [removeme index];
                    end
                end
            end
            if ~isempty(removeme)
                varargin(removeme) = [];
            end
            if ~isempty(varargin)
                obj.vargin = varargin;
            end
            if isempty(obj.catalogs)
                obj.catalogs = getcatalogs;
            end

            if nargin
                % first try with the local plate solver
                [obj.result, filename] = obj.solve(filename, 'solve-field', varargin{:});
                % if fails or not installed, use the web service
                if isempty(obj.result) && ~isempty(filename)
                    obj.solve(filename, 'web', varargin{:});
                end
                % image(obj);
            end

        end

        function obj = local(obj, filename, varargin)
            % LOCAL Loads an image and identifies its objects using local solve-field
            %
            % as = local(astrometry, file, ...);
            %   Solve the given astrophotography image with local method.
            %
            % input(optional):
            %    filename: an image to annotate
            %    any name/value pair as:
            %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
            %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
            %      radius:  approximate field size     (in deg)
            %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
            %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
            %      object:  name of an object on field (string, e.g. 'M 16')
            %
            % Example:
            %   as=local(astrometry, 'M33.jpg','scale-low', 0.5, 'scale-high',2);
            if isempty(varargin) && ~isempty(obj.vargin)
                varargin = obj.vargin;
            end
            [ret, filename] = solve(obj, filename, 'solve-field', varargin{:});
        end

        function obj = web(obj, filename, varargin)
            % WEB Loads an image and identifies its objects using web service
            %
            % as = web(astrometry, file, ...);
            %   Solve the given astrophotography image with web method.
            % web(as)
            %   Once solved, the field is displayed on http://www.sky-map.org
            %
            % input(optional):
            %    filename: an image to annotate
            %    any name/value pair as:
            %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
            %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
            %      radius:  approximate field size     (in deg)
            %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
            %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
            %      object:  name of an object on field (string, e.g. 'M 16')
            %
            % Example:
            %   as=web(astrometry, 'M33.jpg','scale-low', 0.5, 'scale-high',2);
            if nargin == 1 && ~isempty(obj.result) && isstruct(obj.result)
                % display a Sky-Map.org view of the astrometry field
                sz = max([obj.result.RA_max-obj.result.RA_min obj.result.Dec_max-obj.result.Dec_min]);
                z  = 160.*2.^(0:-1:-8); % zoom levels in deg in sky-map
                z  = find(sz*4 > z, 1);
                if isempty(z)
                    z = 9;
                end
                url = sprintf(['http://www.sky-map.org/?ra=%f&de=%f&zoom=%d' ...
                    '&show_grid=1&show_constellation_lines=1' ...
                    '&show_constellation_boundaries=1' ...
                    '&show_const_names=0&show_galaxies=1&img_source=DSS2'], ...
                    obj.result.RA/15, obj.result.Dec, z);
                % open in system browser
                open_system_browser(url);
            else
                if isempty(varargin) && ~isempty(obj.vargin)
                    varargin = obj.vargin;
                end
                if nargin < 2 || isempty(filename)
                    filename = obj.filename;
                end
                [ret, filename] = solve(obj, filename, 'web', varargin{:});
            end
        end

        function [ret, filename] = solve(obj, filename, method, varargin)
            % SOLVE Solve an image field. Plot further results with IMAGE method.
            %
            %  solve(astrometry, filename, method, ...)
            %
            %  input:
            %    filename: an image to annotate
            %    method:   'solve-field' (default) or 'nova'
            %    any name/value pair as:
            %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
            %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
            %      radius:  approximate field size     (in deg)
            %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
            %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
            %      object:  name of an object on field (string, e.g. 'M 16')
            %
            % example: as.solve('M33.jpg')
            %          as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','radius', 2)
            %          as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','radius', 2,...
            %             'scale-low',0.5,'scale-high',2)

            % varargin: ra, dec, radius, scale-low, scale-high

            ret = [];

            if nargin < 3,      method = ''; end
            if isempty(method), method = 'solve-field';
            end

            if ishold(obj)
                disp(' Current solver is RUNNING.')
                return
            end

            isnova = strcmp(method, 'nova') || strcmp(method, 'web') || strcmp(method, 'astrometry.net') || strcmp(method, 'client.py');

            % check if executables are available
            if isnova
                if isempty(obj.executables.client_py) || isempty(obj.executables.python)
                    return
                end
            else
                if isempty(obj.executables.solve_field)
                    return
                end
            end

            % request image if missing
            if nargin < 2,        filename = '';            end
            if isempty(filename), filename = obj.filename; end
            if iscell(filename)
                ret = {};
                for index = 1:numel(filename)
                    ret{end+1} = solve(obj, filename{index}, method, varargin{:});
                end
                return
            end
            if ischar(filename) && isempty(filename)
                % request an image file to solve
                [filename, pathname] = uigetfile( ...
                    {'*.JPG;*.JPEG;*.jpg;*.jpeg;*.GIF;*.gif;*.PNG;*.png;*.FITS;*.fits;*.FTS;*.fts' 'All supported images (JPG,PNG,GIF,FITS)';
                    '*.JPG;*.JPEG;*.jpg;*.jpeg' 'JPEG image (JPG)';
                    '*.GIF;*.gif' 'GIF image (GIF)';
                    '*.PNG;*.png' 'PNG image (PNG)';
                    '*.FITS;*.fits;*.FTS;*.fts' 'FITS image (FITS)';
                    '*.*',  'All Files (*.*)'}, ...
                    ' Pick an astrophotography image to solve');
                if isequal(filename,0)
                    filename = '';
                    return
                end
                filename = fullfile(pathname, filename);
            elseif isempty(filename)
                return % nothing to do
            end
            % when given as RGB
            if isnumeric(filename)
                rgb = filename;
                filename = tempname;
                imwrite(im, [filename '.png'], 'png', 'Author', getenv('USER'), 'CreationTime', datestr(now), 'Software', mfilename);
            end
            % handle distant images
            if ischar(filename) && any(strncmp(filename, {'http' 'ftp:'}, 4))
                [p,f,e] = fileparts(filename);
                tmpfile = [tempname e];
                filename = urlwrite(filename, tmpfile);
            end
            if ~ischar(filename)
                filename = '';
                return
            elseif any(strcmpi(filename,{'/dev/null' 'null' 'none'}))
                filename = [];
                return
            elseif isempty(dir(filename))
                disp([' ERROR: invalid file name "' filename '"'])
                filename = [];
                return
            end
            obj.filename = filename;

            % build the command line
            if isempty(obj.process_dir) % first call
                obj.process_dir = tempname;
            elseif ~isempty(dir(obj.process_dir))
                % clean any previous annotation
                rmdir(obj.process_dir, 's');
            end
            if ~isfold(obj.process_dir)
                mkdir(obj.process_dir);
            end
            d=obj.process_dir;

            if isnova
                % is there an API_KEY ? request it if missing...
                if isempty(obj.api_key) && isempty(getenv('AN_API_KEY'))
                    % request the API_KEY via a dialogue
                    op.Resize      = 'on';
                    op.WindowStyle = 'normal';
                    op.Interpreter = 'tex';
                    prompt = ['{\color{blue}Enter a nova.astrometry.net API key}' 10 ...
                        '  (e.g. "slsratwckfdyxhjq")' 10 ...
                        'Create an account at {\color{red}http://nova.astrometry.net} (free)' 10 ...
                        'Connect to your account and click on the "{\color{blue}API}" tab to get the key' 10 ...
                        'or define the {\color{blue}AN\_API\_KEY} environment variable'];
                    answer = inputdlg( prompt, 'Input Astrometry.net API key', 1, {'slsratwckfdyxhjq'}, op);
                    if isempty(answer)
                        return
                    end
                    obj.api_key = answer{1};
                end

                cmd = [obj.executables.python ' ' obj.executables.client_py ' --wait'];
                if ~isempty(obj.api_key)
                    cmd = [cmd ' --apikey=' obj.api_key];
                end
                cmd = [cmd ' --upload='   filename];
                cmd = [cmd ' --annotate=' fullfile(d, 'results.json')];
                cmd = [cmd ' --newfits='  fullfile(d, 'results.fits')];
                cmd = [cmd ' --kmz='      fullfile(d, 'results.kml')];
                cmd = [cmd ' --corr='     fullfile(d, 'results.corr')];
            else
                cmd = [obj.executables.solve_field];
                cmd = [cmd ' ' filename];
                if ~isempty(obj.executables.sextractor)
                    cmd = [cmd ' --use-sextractor']; % highly improves annotation efficiency
                end
                cmd = [cmd ' --dir '      d];
                cmd = [cmd ' --new-fits ' fullfile(d, 'results.fits')];
                cmd = [cmd ' --rdls '     fullfile(d, 'results.rdls')];
                cmd = [cmd ' --corr '     fullfile(d, 'results.corr') ' --tag-all'];
                if ~isempty(obj.executables.wcs2kml)
                    cmd = [cmd ' --kmz '      fullfile(d, 'results.kml') ' --no-tweak'];
                end
            end
            cmd = [cmd ' --wcs='      fullfile(d, 'results.wcs')];

            % handle arguments: target name/findobj
            if nargin > 3
                obj = []; rmarg = [];
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
                        if isnova
                            if strcmp(varargin{f}, 'scale-low')
                                varargin{f}='scale-lower';
                            elseif strcmp(varargin{f}, 'scale-high')
                                varargin{f}='scale-upper';
                            end
                        else
                            if strcmp(varargin{f}, 'scale-lower')
                                varargin{f} = 'scale-low';
                            elseif strcmp(varargin{f}, 'scale-upper')
                                varargin{f} = 'scale-high';
                            end
                        end
                        cmd = [cmd ' --' varargin{f} '=' num2str(varargin{f+1})];
                    end
                end
            end

            % execute command
            disp(cmd)
            obj.status    = 'running';
            obj.starttime = clock;
            notify(obj, 'busy');
            disp([' ' method ' please wait (may take e.g. few minutes)...'])

            % create the timer for auto update
            if isempty(obj.timer) || ~isa(obj,timer, 'timer') || ~isvalid(obj.timer)
                obj.timer = timer('TimerFcn', @TimerCallback, 'Period', 5.0, 'ExecutionMode', 'fixedDelay', 'UserData', obj, 'Name', mfilename);
            end
            if strcmp(obj.timer.Running, 'off')
                start(obj.timer);
            end

            % launch a Java asynchronous command
            obj.process_java = java.lang.Runtime.getRuntime().exec(cmd);

            % we shall monitor the completion with a timer
            notify(obj, 'annotationStart');
            ret = cmd;
        end

        function ret = load(obj, d, varargin)
            % LOAD Load astrometry files (WCS,FITS) from a directory
            %   LOAD(astrometry, directory);
            %   The directory may contain WCS, CORR, RDLS or JSON, and image.
            %   No solve plate is performed, only data is read.
            %
            %   LOAD(astrometry, image, ...) starts solve-plate annotation, same as SOLVE.
            if nargin < 2
                d = obj.process_dir;
            end
            if ~isfold(d)
                % loading an image ?
                try
                    im = imread(d);
                    solve(obj, d, varargin{:});
                    return
                end
            end

            obj.result = getresult(d, obj);
            if isempty(obj.result)
                obj.status = 'failed';
                obj.process_java = [];
            else
                obj.status = 'success';
                % is the image available ? use one from the directory
                if ~exist(obj.filename, 'file')
                    % search in the result directory
                    d = [dir(fullfile(obj.result.dir, '*.png')); dir(fullfile(obj.result.dir, '*.fits'))];
                    if ~isempty(d)
                        d = d(1);
                        obj.filename = fullfile(obj.result.dir, d.name);
                    end
                end
            end
            ret = obj.result;
        end

        function ret = getstatus(obj)
            % GETSTATUS Return the astrometry status (success, failed)
            ret = obj.status;
        end

        function st = ishold(obj)
            % ISHOLD get the solver status (IDLE, BUSY)
            %   st = ISHOLD(s) returns 1 when the solver is BUSY.
            st = ~isempty(obj.process_java);
        end

        function fig = plot(obj, varargin)
            % PLOT Show the solve-plate image with annotations. Same as image.
            %
            %   as=astrometry(file);
            %   plot(as);
            fig = image(obj, varargin{:});
        end

        function fig = image(obj, mag)
            % IMAGE Show the solve-plate image with annotations
            %   IMAGE(as, mag) limits the objects up to given magnitude.
            %
            %   as=astrometry(file);
            %   image(as);

            if nargin < 2, mag = inf; end

            fig = [];
            if ~ischar(obj.filename) || isempty(obj.filename) || isempty(dir(obj.filename))
                return
            end
            try
                im  = imread(obj.filename);
            catch ME
                getReport(ME)
                disp([mfilename  ': ERROR: can not read image ' obj.filename])
                return
            end
            fig = figure('Name', [' ' obj.filename]);
            image(im);
            clear im;

            ret = obj.result;

            % set title
            [p,f,e] = fileparts(obj.filename);
            if isfield(ret, 'Constellation')
                title([f e ' in ' ret.Constellation]);
            else
                title([f e]);
            end

            % overlay results
            if ~isempty(obj.result) && strcmp(obj.status, 'success') && isfield(obj.result, 'RA_hms')

                hold on
                % central coordinates
                sz = obj.result.size/2;
                h  = plot(sz(1), sz(2), 'r+'); set(h, 'MarkerSize', 16);
                hcmenu = uicontextmenu;
                uimenu(hcmenu, 'Label', '<html><b>Field center</b></html>');
                uimenu(hcmenu, 'Label', ['RA=  ' ret.RA_hms]);
                uimenu(hcmenu, 'Label', ['DEC= ' ret.Dec_dms]);
                uimenu(hcmenu, 'Label', ['Rotation= ' num2str(ret.rotation) ' [deg]']);
                set(h, 'UIContextMenu', hcmenu);

                % get list of visible objects
                v = visible(obj, mag);

                % make sure we do not plot too many
                if numel(v) > 1500
                    v = v(1:1500);
                end

                for index = 1:numel(v)
                    % find all objects from data base within bounds
                    this = v(index);

                    % stars in green, DSO in cyan
                    if strcmp(this.catalog,'stars')
                        c = 'g';
                    else
                        c = 'c';
                    end
                    x = this.X; y = this.Y;

                    % plot symbol
                    if isfinite(this.SIZE) && this.SIZE > 5
                        h = plot(x, y, [c 'o']);
                        set(h, 'MarkerSize', ceil(this.SIZE));
                    else
                        h = plot(x, y, [c 's']); this.SIZE=12;
                    end

                    % context menu
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
                        t=text(x+this.SIZE,y-this.SIZE,this.NAME); set(t,'Color', c);
                    end
                end
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

        function found = findobj(obj, name)
            % Find a given object in catalogs.
            catalog_names = fieldnames(obj.catalogs);
            found = [];

            % check first for name without separator
            if ~any(name == ' ')
                [n1, n2] = strtok(name, '0123456789');
                found = findobj(obj, [n1 ' ' n2]);
                if ~isempty(found)
                    return
                end
            end
            namel = strtrim(lower(name));
            for f = catalog_names(:)'
                catalog = obj.catalogs.(f{1});
                if ~isfield(catalog, 'MAG')
                    continue
                end
                NAME = lower(catalog.NAME);
                NAME = regexprep(NAME, '\s*', ' ');
                % search for name
                index = find(~cellfun(@isempty, strfind(NAME, [';' namel ';'])));
                if isempty(index)
                    index = find(~cellfun(@isempty, strfind(NAME, [namel ';'])));
                end
                if isempty(index)
                    index = find(~cellfun(@isempty, strfind(NAME, [';' namel])));
                end
                if isempty(index)
                    index = find(~cellfun(@isempty, strfind(NAME, namel)));
                end
                if ~isempty(index)
                    found.index   = index(1);
                    found.catalog = f{1};
                    found.RA      = catalog.RA(found.index);
                    found.DEC     = catalog.DEC(found.index);
                    found.MAG     = catalog.MAG(found.index);
                    found.TYPE    = catalog.TYPE{found.index};
                    found.NAME    = catalog.NAME{found.index};
                    found.DIST    = catalog.DIST(found.index);
                    break
                end
            end

            if ~isempty(found)
                disp([mfilename  ': Found object ' name ' as: ' found.NAME])
                if found.DIST > 0
                    fprintf('  %s: Magnitude: %.1f ; Type: %s ; Dist: %.3g [ly]\n', found.catalog, found.MAG, found.TYPE, found.DIST*3.262)
                else
                    fprintf('  %s: Magnitude: %.1f ; Type: %s\n', found.catalog, found.MAG, found.TYPE )
                end
                if ~isempty(obj.result) && isfield(obj.result, 'RA_hms')
                    ret = obj.result;
                    if ret.RA_min <= found.RA && found.RA <= ret.RA_max && ret.Dec_min<= found.DEC && found.DEC <= ret.Dec_max
                        disp(['    ' found.NAME ' is within the image field.'])
                    end
                end
            else
                disp([mfilename  ': object ' name ' was not found.'])
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
                    disp(obj.filename)
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
            fprintf('%s = %s for "%s":\n', iname, id, obj.filename)
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
            elseif isfold(obj.process_dir)
                if isdeployed || ~usejava('jvm') || ~usejava('desktop')
                    disp(['  ' upper(obj.status) ' in ' obj.process_dir]);
                else
                    disp(['  ' upper(obj.status) ' in <a href="' obj.process_dir '">' obj.process_dir '</a>']);
                end
            else
                disp(['  ' upper(obj.status) ': use annotate(as,''filename'').']);
            end

        end

        function display(obj)
            % DISPLAY Display Astrometry object (short)

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
            if ishold(obj)
                fprintf('%s = %s for "%s" BUSY\n',iname, id, obj.filename)
            else
                fprintf('%s = %s for "%s"\n',iname, id, obj.filename)
            end
        end

        function waitfor(obj)
            % WAITFOR Waits for completion of the annotation
            while ishold(obj)
                pause(5)
            end
        end

        function stop(obj)
            % STOP Ends any current annotation and reset the object.
            % clear the timer
            if ~isempty(obj.timer) && isa(obj.timer, 'timer')
                stop(obj.timer);
                delete(obj.timer);
            end
            obj.timer = [];
            p = obj.process_java;
            if ~isempty(p) && isjava(p)
                try
                    p.destroy;
                    disp([mfilename  ': abort current annotation...'])
                end
            end
            obj.process_java = [];
            obj.status = 'failed';
        end

        function c = get_catalogs(obj)
            % GET_CATALOGS Get the loaded catalogs
            c = obj.catalogs;
        end

        function st = get_state(obj)
            % GET_STATE Return the astrometry state, e.g. BUSY, FAILED, SUCCESS.
            st = obj.status;
        end

    end

end


function ret = getresult(d, obj)
% getresult: extract WCS and star matching information from the output files.
%
% input:
%   d: directory where astrometry.net results are stored.

ret = [];
for file = {'results.wcs' 'wcs.fits'}
    if exist(fullfile(d, file{1}), 'file')
        ret.wcs  = read_fits(fullfile(d, file{1}));

        % get image center and print it
        if isfield(ret.wcs,'meta') && isfield(ret.wcs.meta,'CRVAL1')
            wcs = ret.wcs.meta;
            ret.wcs.meta.CD = [wcs.CD1_1 wcs.CD1_2 ; wcs.CD2_1 wcs.CD2_2];

            % get central coordinates
            ret.size = [wcs.IMAGEW wcs.IMAGEH];
            sz  = ret.size/2;

            [ret.RA, ret.Dec] = xy2sky_tan(ret.wcs.meta, sz(1), sz(2)); % MAAT Ofek (private)
            ret.RA       = rad2deg(ret.RA);
            ret.Dec      = rad2deg(ret.Dec);
            ret.RA_hms   = getra(ret.RA/15,true);
            ret.Dec_dms  = getdec(ret.Dec, true);
            
            ret.pixel_scale = sqrt(abs(wcs.CD1_1 * wcs.CD2_2  - wcs.CD1_2 * wcs.CD2_1))*3600; %pixel scale (arcsec/pixel)
            
            ret.rotation = atan2d(wcs.CD2_1, wcs.CD1_1); %rotation angle

            % compute RA,Dec image bounds
            RA = []; Dec = [];
            [RA(end+1), Dec(end+1)] = xy2sky_tan(ret.wcs.meta, wcs.IMAGEW, wcs.IMAGEH);
            [RA(end+1), Dec(end+1)] = xy2sky_tan(ret.wcs.meta, 1         , wcs.IMAGEH);
            [RA(end+1), Dec(end+1)] = xy2sky_tan(ret.wcs.meta, wcs.IMAGEW, 1);
            [RA(end+1), Dec(end+1)] = xy2sky_tan(ret.wcs.meta, 1         , 1);
            RA  = rad2deg(RA);
            Dec = rad2deg(Dec);
            ret.RA_min  = min(RA);
            ret.RA_max  = max(RA);
            ret.Dec_min = min(Dec);
            ret.Dec_max = max(Dec);

            % identify constellation we are in
            [m, index] = min( (ret.RA - obj.catalogs.constellations.RA).^2 + (ret.Dec- obj.catalogs.constellations.DEC).^2 );
            ret.Constellation = obj.catalogs.constellations.Name{index};
        end
    end
end
for file = {'results.rdls' 'rdls.fits'}
    if exist(fullfile(d, file{1}), 'file')
        ret.rdls = read_fits(fullfile(d, file{1}));
    end
end
for file = {'results.corr' 'corr.fits'}
    if exist(fullfile(d, file{1}), 'file')
        ret.corr = read_fits(fullfile(d, file{1}));
    end
end
if exist(fullfile(d, 'results.json'), 'file')
    ret.json = loadjson(fullfile(d, 'results.json'));
end
if ~isempty(ret)
    ret.dir     = d;
end
end


function executables = find_executables
% find_executables: locate executables, return a structure

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
    for exe = ["solve-field" "sextractor" "python" "python3" "wcs2kml" "client.py"] %try these programs
        for cmd = [fullfile(root_fold,exe+ext) fullfile(root_fold,exe) exe+ext exe] %try these calls
            if isfile(cmd)
                err = 0;
            else
                [err,~] = system(cmd_prefix + cmd + " --version"); %run from Matlab
            end
            name = regexprep(exe, {'-' '\.'}, '_');
            if strcmp(name, 'python3')
                name = 'python';
            end
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


function data = getcatalogs
%Load star catalogs
persistent old_data
if ~isempty(old_data)
    data = old_data;
else
    fprintf('Loading star catalogs from %s.mat\n',mfilename) %progress
    data = load(mfilename); %load data
    for f = fieldnames(data)'
        name = f{1};
        if ~isempty(data.(name))
            num  = numel(data.(name).RA);
            if isfield(data.(name), 'Description')
                desc = data.(name).Description;
            else
                desc = '';
            end
            fprintf(' %s (%g entries) - %s\n',name,num,desc) %display info
        end
    end
    old_data = data; %cache for next time
end
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


function TimerCallback(src, ~)
% TimerCallback: update status/view from timer event
obj = get(src, 'UserData');
if isvalid(obj)
    try
        % check if any astrometry job is running
        exitValue = 0;
        if ~isempty(obj.process_java) && isjava(obj.process_java)
            try
                exitValue = obj.process_java.exitValue; % will raise error if process still runs
                active  = 0;
            catch ex %#ok<NASGU>
                % still running
                if isempty(obj.process_java) || ~isjava(obj.process_java)
                    active  = 0;
                else
                    active  = 1;
                end
            end
            % not active anymore: process has ended.
            if ~active
                if exitValue ~= 0
                    obj.result = [];
                    obj.status = 'failed';
                else
                    load(obj)
                    if ~isempty(obj.result) obj.status = 'success';
                    else obj.status = 'failed';
                    end
                end
                disp([' annotation end: ' upper(obj.status) ' for ' obj.filename '. exit value=' num2str(exitValue)]);
                % clear the timer
                stop(src); delete(src); obj.timer = [];
                obj.process_java = [];
                obj.duration = etime(clock, obj.starttime);

                notify(obj, 'annotationEnd');
                notify(obj, 'idle');
                if obj.autoplot
                    image(obj);
                    assignin('base', 'ans', obj);
                    out = obj;
                    disp(out)
                end
                beep
            end
        end
    catch ME
        getReport(ME)
    end
else
    delete(src);
end
end