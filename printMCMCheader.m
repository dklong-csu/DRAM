function printMCMCheader(fileID)
%printMCMCheader Print output header
%   TODO write documentation

%   Handle defaults arguments
arguments
    %   fileID=1 prints to the command window
    %   fileID=2 creates an error message
    %   Other scalars print to a file associated with that identifier
    fileID {mustBeScalarOrEmpty} = 1
end

%   Horizontal line on top
%       Step
fprintf(fileID,"\x250F");
for vvv=1:10
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x2533");

%       Acceptance rate
for vvv=1:17
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x2533");

%       Rhat
for vvv=1:12
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x2533");

%       Effective sample size
for vvv=1:10
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x2533");

%       Time since last update
for vvv=1:17
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x2513");
fprintf(fileID,"\n");

%   Column headers

fprintf(fileID,"\x2503");
fprintf(fileID,"%-10s"," Step");
fprintf(fileID,"\x2503");

fprintf(fileID,"%-16s "," Acceptance Rate");
fprintf(fileID,"\x2503");

fprintf(fileID,"%-12s"," Worst Rhat");
fprintf(fileID,"\x2503");

fprintf(fileID,"%-10s","   ESS");
fprintf(fileID,"\x2503");

fprintf(fileID,"%-17s"," Time (mm:ss.SSS)");
fprintf(fileID,"\x2503");

fprintf(fileID,"\n");

%   Extra information for user
%   Step - no recommendation
fprintf(fileID,"\x2503");
fprintf(fileID,"%10s"," ");
fprintf(fileID,"\x2503");

%   Acceptance rate: show median and 95% interval
fprintf(fileID,"%17s"," 2.5%-50%-97.5% ");
fprintf(fileID,"\x2503");

%   Rhat: recommend less than 1.01
fprintf(fileID,"%12s","  <1.01    ");
fprintf(fileID,"\x2503");

%   ESS: minimum 50 per chain
fprintf(fileID,"%10s"," >50    ");
fprintf(fileID,"\x2503");

%   Time - n/a
fprintf(fileID,"%17s"," ");
fprintf(fileID,"\x2503");

fprintf(fileID,"\n");

%   Horizontal line
%       Step
fprintf(fileID,"\x2503");
for vvv=1:10
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x254B");

%       Acceptance rate
for vvv=1:17
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x254B");

%       Rhat
for vvv=1:12
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x252B");

%       Effective sample size
for vvv=1:10
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x252B");

%       Time since last update
for vvv=1:17
    fprintf(fileID,"\x2501");
end
fprintf(fileID,"\x252B");

fprintf(fileID,"\n");
end