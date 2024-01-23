function printMCMCupdate(step,n_accepted,rhatvec,ESS,elapsed_time,fileID)
%printMCMCupdate Prints an update of the MCMC iterations
%   TODO document

%   Handle defaults arguments
arguments
    step {mustBeInteger} = 1
    n_accepted {mustBeInteger} = 0
    rhatvec {mustBeFloat} = Inf
    ESS {mustBeFloat} = 0
    %   FIXME Not sure how to do this one
    elapsed_time = -(datetime('now') - datetime('now'));
    %   fileID=1 prints to the command window
    %   fileID=2 creates an error message
    %   Other scalars print to a file associated with that identifier
    fileID {mustBeScalarOrEmpty} = 1
end

%   Step
fprintf(fileID,"\x2503");
fprintf(fileID," %-9d",step);
fprintf(fileID,"\x2503");

%   Acceptance rate
AR = 100*n_accepted/step;
q = quantile(AR,[0.025, 0.5, 0.975]);
fprintf(fileID," %3.0f%% %3.0f%% %3.0f%%  ",q);
fprintf(fileID,"\x2503");

%   Rhat
fprintf(fileID,"% 9.3f   ",max(rhatvec));
fprintf(fileID,"\x2503");

%   Effective sample size
fprintf(fileID,"% -10d",max(0,ceil(min([ESS,step]))));
fprintf(fileID,"\x2503");

%   Elapsed time
elapsed_time.Format = 'mm:ss.SSS';
fprintf(fileID,"%17s",elapsed_time);
fprintf(fileID,"\x2503");

fprintf(fileID,"\n");
end