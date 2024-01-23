clc
close all
clear variables
delete(gcp("nocreate"))

rng(1)
%%  Test DRAM

%   Settings
x0 = mvnrnd([0,0,0],10*eye(3),6)';


%   Run algorithm
mcmcData = DRAM(@rosenbrock, @prior, x0,...
    "delayStages",5,...
    "nSamples",10000,...
    "nBurnedSamples",100,...
    "nAdaptSamples",100,...
    "updateEveryN",100,...
    "nChains",6);

%%   Generate plot
mymap = [1,1,1;parula];
vars = ["Var1","Var2","Var3"];
t = plotPosterior(mcmcData,vars,mymap);






%%  Test functions
function y = rosenbrock(x)
    d = length(x);
    y = mvnpdf(resize(x,[1,d]),zeros(1,d),eye(d));
    y = log(y);
end

function y = prior(x)
    y = 0;
end
