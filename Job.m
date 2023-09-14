%Class for starting and monitorying system commands.
%
%Example:
% clc, J=Job, pause(1); J.Run('sleep 10'), J, pause(1); J

classdef Job < handle

    properties (SetAccess = protected) %read only
        status   = '' %can be: init, running, failed, success
        start    = []
        duration = duration;
        worker   = []
        timer    = []
    end

    events
        idle
        starting
        busy
        stopping
    end

    methods

        function obj = Job(cmd) % Class construtor
            if nargin && ~isempty(cmd)
                obj.Run(cmd)
            end
        end

        function Run(obj,cmd) % Execute command
            obj.start = datetime;
            obj.status = 'running';
            obj.notify('busy');
            if isempty(obj.timer)
                obj.timer = timer(TimerFcn=@(o,e)CheckStatus(o,e,obj), ExecutionMode='fixedDelay', Period=0.1); %timer for auto update
            end
            if strcmp(obj.timer.Running, 'off')
                obj.timer.start
            end
            obj.worker = java.lang.Runtime.getRuntime().exec(cmd); %start job
            obj.notify('starting'); % Wait for completion
        end

        function out = getstatus(obj) %status: 'success' 'failed'
            out = obj.status;
        end

        function Wait(obj) % Waits till job is finished
            while ~isempty(obj.worker)
                pause(1)
            end
        end

        function Stop(obj) % Abort the Job
            if ~isempty(obj.timer)
                stop(obj.timer);
                delete(obj.timer);
            end
            obj.timer = [];
            p = obj.worker;
            if ~isempty(p) && isjava(p)
                try
                    p.destroy;
                end
            end
            obj.worker = [];
            obj.status = 'failed';
        end
        function st = get_state(obj) %states: BUSY, FAILED, SUCCESS
            st = obj.status;
        end
    end
end

function CheckStatus(src, x,obj) % Update job status
% obj = get(src, 'UserData');
exitValue = 0;
if ~isempty(obj.worker) && isjava(obj.worker)
    obj.duration = datetime - obj.start;
    try
        obj.worker
        exitValue = obj.worker.exitValue; %error if process still runing
        active = 0;
    catch ex %#ok<NASGU>
        % still running
        if isempty(obj.worker) || ~isjava(obj.worker)
            active = 0;
        else
            active = 1;
        end
    end
    
    % Close if process has ended
    if ~active
        if exitValue ~= 0
            % obj.result = [];
            obj.status = 'failed';
        else
        %     if ~isempty(obj.result)
        %         obj.status = 'success';
        %     else
        %         obj.status = 'failed';
        %     end
        end
        
        % Clean up
        src.stop
        src.delete
        obj.timer = [];
        obj.worker = [];
        obj.notify('stopping');
        obj.notify('idle');

        disp('Done')
    end
end
end