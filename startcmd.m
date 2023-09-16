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