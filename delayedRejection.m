function [yn, accept, ynpost, endStage] = delayedRejection(xc, Covar, logPosterior, delayChance, logPostxc, maxStage)
    %   Set up for delayed rejections
    stage = 1;
    [yn, accept, ynpost, endStage] = recursiveDR(xc, xc, Covar, logPosterior, delayChance, -Inf, logPostxc, stage, maxStage);
end


%%  Helper function
function [yn, accept, ynpost, endStage] = recursiveDR(xc, y, Covar, logPosterior, delayChance, postybest, lpxc, stage, maxStage)
    %   Draw from distribution
    yi = mvnrnd(y,Covar,1);
    %   Evaluate at proposed point
    lp = logPosterior(yi);
    %   Compute acceptance logic
    Ni = log( max(0, exp(lp) - exp(postybest) ) );
    Di = log( exp(lpxc) - exp(postybest) );
    ratio = Ni - Di;
    r = log( unifrnd(0,1) );
    accept = (ratio >= r);
    %   Either accept or go through delayed rejection logic
    if accept
        yn = yi;
        accept = true;
        ynpost = lp;
        endStage = stage;
    else
        stage = stage+1;
        %   End because probability of not delaying
        r = unifrnd(0,1);
        if r > delayChance || stage > maxStage
            %   Then we do not allow a delayed rejection
            accept = false;
            yn = xc;
            ynpost = lp;
            endStage = stage;
        else
            %   Recursion
            scale = 0.5;
            postybest = max(postybest, lp);
            [yn, accept,ynpost, endStage] = recursiveDR(xc, xc, scale*Covar, logPosterior, delayChance, postybest, lpxc, stage, maxStage);
        end
    end

end