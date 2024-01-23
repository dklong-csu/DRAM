function [rhat,B,W,varp,svar] = calcRhat(samples)
    chain_avg = mean(samples,1);
    
    [N,M] = size(samples);

    %   Between chain variance
    B = N * var(chain_avg);

    %   Within chain variance
    svar = var(samples);
    W = mean( svar );

    %   Var+ estimator
    varp = (N-1)/N * W + B/N;

    %   Rhat ratio
    rhat = sqrt( varp/W );
end



