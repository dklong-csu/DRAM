function mcmcData = DRAM(logLikelihood, logPrior, x0, options)
%DRAM Delayed Rejection Adaptive Metropolis
%   TODO documentation

%%   Initialization
arguments
    %   FIXME
    logLikelihood {mustBeNonempty}
    logPrior {mustBeNonempty}
    x0 {mustBeNumeric}

    %   Number of samples to generate
    options.nBurnedSamples = 50
    options.nAdaptSamples = 100
    options.nSamples = 200
    options.nChains = size(x0,2)

    options.nDims {mustBeNonnegative,mustBeInteger} = size(x0,1)
    options.useParallel = true;

    %   DRAM settings
    options.initCov {mustBeNumeric} = eye(size(x0,1))
    options.covarReduce {mustBePositive,mustBeFloat} = 0.5
    options.delayStages {mustBeInteger,mustBeNonnegative} = 2
    options.covEpsilon {mustBeNonnegative, mustBeFloat} = 0

    %   Display
    options.verbose {mustBeNumericOrLogical} = true;
    options.fileID = 1;
    options.updateEveryN {mustBeInteger,mustBeNonnegative} = 10;
    options.reprintHeaderEveryN {mustBeInteger,mustBeNonnegative} = 25;

    %   Output data settings
    options.dataFilename {mustBeTextScalar} = "myDRAMdata.mat";
end

if options.useParallel
    p = parpool('Threads');
end

%   Data object to store results
mcmcData.options = options;

%   Set up chains
%   Matlab stores data in column-major format which means incrementing the
%   first index moves incrementally in memory. The most demanding
%   computations on the samples will be on vectors representing a single
%   parameter dimension over many iterations. Thus we store the samples in
%   an array following
%       samples( iteration, parameter dimension, chain)
%   Then, extracting the vector samples(:,p,c) will give the entire set of
%   samples for a single parameter for a single chain in contiguous memory.
%   Similarly, extracting the matrix samples(:,:,c) will give the entire set
%   of samples of a single chain in contiguous memory.

burnedSamples = zeros(options.nBurnedSamples+1, options.nDims, options.nChains);
burnedLogPost = zeros(options.nBurnedSamples+1, options.nChains);
burnedAccepts = zeros(options.nChains,1);

adaptSamples = zeros(options.nAdaptSamples+1, options.nDims, options.nChains);
adaptLogPost = zeros(options.nAdaptSamples+1, options.nChains);
adaptAccepts = zeros(options.nChains,1);

samples = zeros(options.nSamples+1, options.nDims, options.nChains);
samplesLogPost = zeros(options.nSamples+1, options.nChains);
nAcceptances = zeros(options.nChains,1);

%   Initializes the samples data, but this will need to be updated once the
%   samples are actually generated
mcmcData.burnedSamples = burnedSamples;
mcmcData.burnedLogPost = burnedLogPost;
mcmcData.burnedAccepts = burnedAccepts;

mcmcData.adaptSamples = adaptSamples;
mcmcData.adaptLogPost = adaptLogPost;
mcmcData.adaptAccepts = adaptAccepts;

mcmcData.samples = samples;
mcmcData.samplesLogPost = samplesLogPost;
mcmcData.nAcceptances = nAcceptances;

mcmcData.ESS = 0;
mcmcData.IAT = 1;
mcmcData.rhat = Inf*ones(options.nDims,1);

logPosterior = @(x) logLikelihood(x) + logPrior(x);
covarianceReductionFactor = options.covarReduce;
numDelayStages = options.delayStages;

%   Initial evaluation for each chain
parfor ccc=1:options.nChains
    %   TODO error checking to make sure enough x0
    lp = logPosterior(x0(:,ccc));
    burnedSamples(1,:,ccc) = x0(:,ccc);
    burnedLogPost(1,ccc) = lp;
end

%%  Burn in phase
%   Burn some samples to find a better starting point (if necessary)
if options.verbose && options.nBurnedSamples>0
    fprintf(options.fileID, "\n%-s\n","Burned Samples Summary");
    printMCMCheader(options.fileID);
end

burnedCov = options.initCov;
for iii=2:options.nBurnedSamples+1
    %   Time full iteration
    if iii==2
        t1 = datetime('now');
    end

    %   Do each chain independently
    xc = squeeze(burnedSamples(iii-1,:,:));
    xcPost = burnedLogPost(iii-1,:);
    parfor ccc=1:options.nChains
        %   Delayed Rejection sampling
        [y, accept, ypost, ~] = delayedRejection( xc(:,ccc), ...
            burnedCov, ...
            logPosterior,...
            covarianceReductionFactor,...
            xcPost(ccc),...
            numDelayStages);
        %   Update arrays
        burnedSamples(iii, :, ccc) = y;
        burnedLogPost(iii,ccc) = ypost;
        burnedAccepts(ccc) = burnedAccepts(ccc) + accept;
    end

    %   Update every N
    if options.verbose && or(rem(iii-1,options.updateEveryN)==0,iii-1==options.nBurnedSamples)
        %   Occasionally reprint the header so you don't have to scroll up
        %   far to remember what each column is
        if rem( (iii-1)/options.updateEveryN, options.reprintHeaderEveryN )==0
            printMCMCheader(options.fileID);
        end

        %   FIXME -- this is to fix issues with computing statistics with
        %   few samples
        %   FIXME -- make this a function
        if iii<5
            rhat = Inf*ones(1,options.nDims);
            ESS = zeros(1,options.nDims);
        else
            %   Loop over each dim
            for ddd=1:options.nDims
                %   Split each chain in half
                %   If odd number of samples, ignore very first sample
                if rem(iii,2) == 0
                    mysamps = squeeze(burnedSamples(1:iii,ddd,:));
                else
                    mysamps = squeeze(burnedSamples(2:iii,ddd,:));
                end
                half = size(mysamps,1)/2;
                splitSamps = zeros(half, 2*size(mysamps,2));
                for ccc=1:options.nChains
                    idx1 = 2*ccc - 1;
                    idx2 = 2*ccc;
                    splitSamps(:,idx1:idx2) = reshape(mysamps(:,ccc), [half, 2]);
                end

                %   Compute Rhat
                [rhat(ddd), B, W, varp, svar] = calcRhat(splitSamps);

                %   Effective sample size
                acf = zeros(half, 2*options.nChains);
                for ccc=1:size(splitSamps,2)
                    acf(:,ccc) = autocorr(splitSamps(:,ccc), "NumLags",half-1);
                end
                rho = 1 - (W - mean(svar.*acf,2))/varp;
                P = rho(1:2:end-1) + rho(2:2:end);
                tk = find(P<=0,1);
                if isempty(tk)
                    tau = -1 + 2*sum(P);
                else
                    tau = -1 + 2*sum( P(1:tk-1) );
                end
                ESS(ddd) = numel(splitSamps)/tau;
            end
        end
        %   Print summary
        t2 = datetime('now');
        elapsed_time = t2 - t1;
        printMCMCupdate(iii-1, burnedAccepts, rhat, ESS, elapsed_time, options.fileID);

        %   Update and save current
        mcmcData.burnedSamples = burnedSamples;
        mcmcData.burnedLogPost = burnedLogPost;
        mcmcData.burnedAccepts = burnedAccepts;
        mcmcData.ESS = min(ESS);
        mcmcData.IAT = tau;
        mcmcData.rhat = rhat;
        save(options.dataFilename,"mcmcData");

        %   Reset timing
        t1 = datetime('now');
    end

   
end

%   Finalize this part of the algorithm
if options.verbose
    printMCMCcloser(options.fileID);
end
mcmcData.burnedSamples = burnedSamples;
mcmcData.burnedLogPost = burnedLogPost;
mcmcData.burnedAccepts = burnedAccepts;
save(options.dataFilename,"mcmcData");


%%  Covariance adaptation phase
%   Generate samples prior to adapting the covariance matrix to get a
%   better estimate
if options.verbose && options.nAdaptSamples>0
    fprintf(options.fileID, "\n%-s\n","Adaptation Samples Summary");
    printMCMCheader(options.fileID);
end

%   Use the same initial covariance matrix
adaptationCov = options.initCov;
%   Start from the final sample in the burn in phase
adaptSamples(1,:,:) = mcmcData.burnedSamples(end,:,:);
adaptLogPost(1,:) = mcmcData.burnedLogPost(end,:);

for iii=2:options.nAdaptSamples+1
    %   Time full iteration
    if iii==2
        t1 = datetime('now');
    end

    %   Do each chain independently
    xc = squeeze(adaptSamples(iii-1,:,:));
    xcPost = adaptLogPost(iii-1,:);
    parfor ccc=1:options.nChains
        %   Delayed Rejection sampling
        %   FIXME on covariance reduction factor
        [y, accept, ypost, ~] = delayedRejection( xc(:,ccc), ...
            adaptationCov, ...
            logPosterior,...
            covarianceReductionFactor,...
            xcPost(ccc),...
            numDelayStages);
        %   Update arrays
        adaptSamples(iii, :, ccc) = y;
        adaptLogPost(iii,ccc) = ypost;
        adaptAccepts(ccc) = adaptAccepts(ccc) + accept;
    end

    %   Update every N
    if options.verbose && or(rem(iii-1,options.updateEveryN)==0,iii-1==options.nAdaptSamples)
        %   Occasionally reprint the header so you don't have to scroll up
        %   far to remember what each column is
        if rem( (iii-1)/options.updateEveryN, options.reprintHeaderEveryN )==0
            printMCMCheader(options.fileID);
        end

        %   FIXME -- this is to fix issues with computing statistics with
        %   few samples
        %   FIXME -- make this a function
        if iii<5
            rhat = Inf*ones(1,options.nDims);
            ESS = zeros(1,options.nDims);
        else
            %   Loop over each dim
            for ddd=1:options.nDims
                %   Split each chain in half
                %   If odd number of samples, ignore very first sample
                if rem(iii,2) == 0
                    mysamps = squeeze(adaptSamples(1:iii,ddd,:));
                else
                    mysamps = squeeze(adaptSamples(2:iii,ddd,:));
                end
                half = size(mysamps,1)/2;
                splitSamps = zeros(half, 2*size(mysamps,2));
                for ccc=1:options.nChains
                    idx1 = 2*ccc - 1;
                    idx2 = 2*ccc;
                    splitSamps(:,idx1:idx2) = reshape(mysamps(:,ccc), [half, 2]);
                end

                %   Compute Rhat
                [rhat(ddd), B, W, varp, svar] = calcRhat(splitSamps);

                %   Effective sample size
                acf = zeros(half, 2*options.nChains);
                for ccc=1:size(splitSamps,2)
                    acf(:,ccc) = autocorr(splitSamps(:,ccc), "NumLags",half-1);
                end
                rho = 1 - (W - mean(svar.*acf,2))/varp;
                P = rho(1:2:end-1) + rho(2:2:end);
                tk = find(P<=0,1);
                if isempty(tk)
                    tau = -1 + 2*sum(P);
                else
                    tau = -1 + 2*sum( P(1:tk-1) );
                end
                ESS(ddd) = numel(splitSamps)/tau;
            end
        end
        %   Print summary
        t2 = datetime('now');
        elapsed_time = t2 - t1;
        printMCMCupdate(iii-1, adaptAccepts, rhat, ESS, elapsed_time, options.fileID);

        %   Update and save current
        mcmcData.adaptSamples = adaptSamples;
        mcmcData.adaptLogPost = adaptLogPost;
        mcmcData.adaptAccepts = adaptAccepts;
        mcmcData.ESS = min(ESS);
        mcmcData.IAT = tau;
        mcmcData.rhat = rhat;
        save(options.dataFilename,"mcmcData");

        %   Reset timing
        t1 = datetime('now');
    end

   
end

%   Finalize this part of the algorithm
if options.verbose
    printMCMCcloser(options.fileID);
end
mcmcData.adaptSamples = adaptSamples;
mcmcData.adaptLogPost = adaptLogPost;
mcmcData.adaptAccepts = adaptAccepts;
save(options.dataFilename,"mcmcData");

%%  Primary sampling phase
%   Burn some samples to find a better starting point (if necessary)
if options.verbose && options.nSamples>0
    fprintf(options.fileID, "\n%-s\n","Primary Samples Summary");
    printMCMCheader(options.fileID);
end

%   Initial covariance matrix computed from the adaptation samples plus a
%   positive perturbation on the diagonal to help ensure a nonsingular
%   matrix
S = zeros((options.nSamples+1)*options.nChains, options.nDims);
for ccc=1:options.nChains
    rowidx1 = (options.nAdaptSamples+1)*(ccc-1)+1;
    rowidx2 = (options.nAdaptSamples+1)*ccc;
    S(rowidx1:rowidx2,:) = mcmcData.adaptSamples(:,:,ccc);
end

eps = options.covEpsilon;
if eps==0
    %   Approximate a decent epsilon based on the determinant of the sample
    %   covariance matrix
    eps = nthroot(1e-6 * det(cov(S)),options.nDims);
end

Cov = adaptCov(options.nDims,S,eps); %fixme
chainCov = zeros(options.nDims, options.nDims, options.nChains);
for ccc=1:options.nChains
    chainCov(:,:,ccc) = Cov;
end

%   Start from the final sample in the burn in phase
samples(1,:,:) = mcmcData.adaptSamples(end,:,:);
samplesLogPost(1,:) = mcmcData.adaptLogPost(end,:);

for iii=2:options.nSamples+1
    %   Time full iteration
    if iii==2
        t1 = datetime('now');
    end

    %   Check if the covariance should be updated now
    if rem(iii-1,options.updateEveryN)==0
        %   Covariance is updated based on all adaptation samples and the
        %   samples in each chain individually
        for ccc = 1:options.nChains
            chainCov(:,:,ccc) = adaptCov(options.nDims, [S;samples(:,:,ccc)], eps);
        end
    end

    %   Do each chain independently
    xc = squeeze(samples(iii-1,:,:));
    xcPost = samplesLogPost(iii-1,:);
    parfor ccc=1:options.nChains
        %   Delayed Rejection sampling
        %   FIXME on covariance reduction factor
        [y, accept, ypost, ~] = delayedRejection( xc(:,ccc), ...
            chainCov(:,:,ccc), ...
            logPosterior,...
            covarianceReductionFactor,...
            xcPost(ccc),...
            numDelayStages);
        %   Update arrays
        samples(iii, :, ccc) = y;
        samplesLogPost(iii,ccc) = ypost;
        nAcceptances(ccc) = nAcceptances(ccc) + accept;
    end

    %   Update every N
    if options.verbose && or(rem(iii-1,options.updateEveryN)==0,iii-1==options.nSamples)
        %   Occasionally reprint the header so you don't have to scroll up
        %   far to remember what each column is
        if rem( (iii-1)/options.updateEveryN, options.reprintHeaderEveryN )==0
            printMCMCheader(options.fileID);
        end

        %   FIXME -- this is to fix issues with computing statistics with
        %   few samples
        %   FIXME -- make this a function
        if iii<5
            rhat = Inf*ones(1,options.nDims);
            ESS = zeros(1,options.nDims);
        else
            %   Loop over each dim
            for ddd=1:options.nDims
                %   Split each chain in half
                %   If odd number of samples, ignore very first sample
                if rem(iii,2) == 0
                    mysamps = squeeze(samples(1:iii,ddd,:));
                else
                    mysamps = squeeze(samples(2:iii,ddd,:));
                end
                half = size(mysamps,1)/2;
                splitSamps = zeros(half, 2*size(mysamps,2));
                for ccc=1:options.nChains
                    idx1 = 2*ccc - 1;
                    idx2 = 2*ccc;
                    splitSamps(:,idx1:idx2) = reshape(mysamps(:,ccc), [half, 2]);
                end

                %   Compute Rhat
                [rhat(ddd), B, W, varp, svar] = calcRhat(splitSamps);

                %   Effective sample size
                acf = zeros(half, 2*options.nChains);
                for ccc=1:size(splitSamps,2)
                    acf(:,ccc) = autocorr(splitSamps(:,ccc), "NumLags",half-1);
                end
                rho = 1 - (W - mean(svar.*acf,2))/varp;
                P = rho(1:2:end-1) + rho(2:2:end);
                tk = find(P<=0,1);
                if isempty(tk)
                    tau = -1 + 2*sum(P);
                else
                    tau = -1 + 2*sum( P(1:tk-1) );
                end
                ESS(ddd) = numel(splitSamps)/tau;
            end
        end
        %   Print summary
        t2 = datetime('now');
        elapsed_time = t2 - t1;
        printMCMCupdate(iii-1, nAcceptances, rhat, ESS, elapsed_time, options.fileID);

        %   Update and save current
        mcmcData.samples = samples;
        mcmcData.samplesLogPost = samplesLogPost;
        mcmcData.nAcceptances = nAcceptances;
        mcmcData.ESS = min(ESS);
        mcmcData.IAT = tau;
        mcmcData.rhat = rhat;
        save(options.dataFilename,"mcmcData");

        %   Reset timing
        t1 = datetime('now');
    end

   
end

%   Finalize this part of the algorithm
if options.verbose
    printMCMCcloser(options.fileID);
end

mcmcData.samples = samples;
mcmcData.samplesLogPost = samplesLogPost;
mcmcData.nAcceptances = nAcceptances;
save(options.dataFilename,"mcmcData");

%%  Final touches
%   Delete the parallel pool
delete(p);
end