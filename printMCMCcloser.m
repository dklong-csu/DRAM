function printMCMCcloser(fileID)
%printMCMCcloser Prints the end of the table iteration summary 
%   TODO document

%   Handle defaults arguments
arguments
    %   fileID=1 prints to the command window
    %   fileID=2 creates an error message
    %   Other scalars print to a file associated with that identifier
    fileID {mustBeScalarOrEmpty} = 1
end

%   Step
fprintf(fileID,"\x2517");
for vvv=1:10
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x253B");

%   Acceptance rate
for vvv=1:17
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x253B");

%   Rhat
for vvv=1:12
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x253B");

%   Effective sample size
for vvv=1:10
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x253B");

%   Run time since last update
for vvv=1:17
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x251B");

fprintf(fileID,"\n");
end