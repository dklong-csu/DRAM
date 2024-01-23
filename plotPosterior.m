function t = plotPosterior(sampsData,varNames,myColorMap)
% plotPosterior - Plots a panel of the posterior distribution from MCMC
% samples.
%
%   TODO: WRITE DOCUMENTATION

%   Extract samples
sampsPre = sampsData.samples;
[N,D,C] = size(sampsPre);
tau = sampsData.IAT;
thin = 1:ceil(tau):N;
M = length(thin);
samps = zeros(M*C,D);
allSamps = zeros(N*C,D);

for ccc=1:C
    r1 = M*(ccc-1)+1;
    r2 = M*ccc;
    samps(r1:r2,:) = sampsPre(thin,:,ccc);
    R1 = N*(ccc-1)+1;
    R2 = N*ccc;
    allSamps(R1:R2,:) = sampsPre(:,:,ccc);
end
%   Extract log posterior
lpPre = sampsData.samplesLogPost;
logPostVals = zeros(M*C,1);
allLogPostVals = zeros(N*C,1);
for ccc=1:C
    r1 = M*(ccc-1)+1;
    r2 = M*ccc;
    logPostVals(r1:r2) = lpPre(thin,ccc);
    R1 = N*(ccc-1)+1;
    R2 = N*ccc;
    allLogPostVals(R1:R2) = lpPre(:,ccc);
end


%   Dimension of parameters
dim = size(samps,2)-1;

%   Some processing to determine axis bounds
buffer = 0.01;
maxprm = max(samps(:,1:end));
minprm = min(samps(:,1:end));
[~,idx] = max(logPostVals);
quants = quantile(samps(:,1:end),[0.025, 0.975]);
n_ticks = 4;
tickformat = "%.1f";


figure
t = tiledlayout(dim,dim,"TileSpacing","compact","Padding","compact");
colormap(myColorMap)

%   Diagonal gives marginal
%   Upper triangle gives scatter plot
%   Lower triangle gives contour plot
for ii=1:dim
    for jj=1:dim
        nexttile(t)
        %   Decide what kind of plot to use
        if ii==jj
            %   Marginal distribution histogram
            h = histogram(samps(:,ii),...
                    'Orientation','vertical');
            hold on
                %   Lines for 95% bounds on marginal
                yplot = linspace(0,max(h.Values));
                xplot1 = quants(1,ii)*(1+0*yplot);
                plot(xplot1,yplot,...
                    '--k',...
                    'LineWidth',2)
                xplot2 = quants(2,ii)*(1+0*yplot);
                plot(xplot2,yplot,...
                    '--k',...
                    'LineWidth',2)
                xplot3 = samps(idx,ii)*(1+0*yplot);
                plot(xplot3,yplot,...
                    '-',...
                    'LineWidth',2,...
                    'Color',[0.8500 0.3250 0.0980]	)
            hold off

            xbuffer = buffer*max( abs( [ minprm(jj),maxprm(jj) ] ) );
            xlim([minprm(ii)-xbuffer, maxprm(ii)+xbuffer])
            tickvec = linspace(minprm(ii), maxprm(ii),n_ticks+2);
            xticks(tickvec(2:end-1))
            xtickformat(tickformat)

            mag = floor( log10( abs( maxprm(ii) ) ) );
            ax = gca;
            if mag <= 1
                ax.XAxis.Exponent = 0;
            else
                ax.XAxis.Exponent = mag;
            end


            ylim([0 max(h.Values)])

            box on
            if ii<dim
                xticklabels('')
            end
            yticklabels('')
            if ii==1
                ylabel(varNames(ii))
            end
            if jj==dim
                xlabel(varNames(jj))
            end
            yticks([])
            set(gca,'LineWidth',2)
        elseif ii>jj
            %   row > column -> lower triangle
            y = allSamps(:,ii);
            x = allSamps(:,jj);
            z = exp(allLogPostVals);
            %   Remove duplicates for the interpolation function
            [uniqueSamples, ia, ~] = unique([x,y],'rows');
            uniqueZ = z(ia);
            margpost = scatteredInterpolant(uniqueSamples(:,1),uniqueSamples(:,2),uniqueZ,'natural','nearest');
        
            xbuffer = buffer*max( abs( [ minprm(jj),maxprm(jj) ] ) );
            ybuffer = buffer*max( abs( [ minprm(ii),maxprm(ii) ] ) );
            plot_range = [minprm(jj) + xbuffer, ...
                maxprm(jj) - xbuffer, ...
                minprm(ii) + ybuffer, ...
                maxprm(ii) - ybuffer];

            fcontour(@(pt1,pt2) margpost(pt1,pt2),...
                plot_range,...
                'Fill','on')
            hold on
                %   95% for x-axis
                xplot1 = quants(1,jj)*ones(1,100);
                xplot2 = quants(2,jj)*ones(1,100);
                yplot1 = linspace(minprm(ii)-ybuffer,maxprm(ii)+ybuffer);
                plot(xplot1,yplot1,...
                    '--k','LineWidth',2)
                plot(xplot2,yplot1,...
                    '--k','LineWidth',2)
                %   95% for y-axis
                yplot2 = quants(1,ii)*ones(1,100);
                yplot3 = quants(2,ii)*ones(1,100);
                xplot3 = linspace(minprm(jj)-xbuffer, maxprm(jj)+xbuffer);
                plot(xplot3, yplot2, ...
                    '--k','LineWidth',2)
                plot(xplot3, yplot3, ...
                    '--k','LineWidth',2)
                %   Show location of best sample
                xbest = samps(idx,jj);
                ybest = samps(idx,ii);
                plot(xbest*ones(1,100),yplot1,...
                    "Color",[0.8500 0.3250 0.0980],...
                    "LineWidth",2)
                plot(xplot3, ybest*ones(1,100),...
                    "Color",[0.8500 0.3250 0.0980],...
                    'LineWidth',2)
                scatter(xbest,ybest,50,[0.8500 0.3250 0.0980],...
                    'filled','o',...
                    'MarkerEdgeColor','k')

            hold off
            box on
            if jj==1
                ylabel(varNames(ii))
            end
            if ii==dim
                xlabel(varNames(jj))
            end
            xlim([minprm(jj)-xbuffer, maxprm(jj)+xbuffer])
            tickvec = linspace(minprm(jj), maxprm(jj), n_ticks+2);
            xticks(tickvec(2:end-1))
            xtickformat(tickformat)

            mag = floor( log10( abs( maxprm(jj) ) ) );
            ax = gca;
            if mag <= 1
                ax.XAxis.Exponent = 0;
            else
                ax.XAxis.Exponent = mag;
            end


            ylim([minprm(ii)-ybuffer, maxprm(ii)+ybuffer])
            tickvec = linspace(minprm(ii), maxprm(ii), n_ticks+2);
            yticks(tickvec(2:end-1))
            ytickformat(tickformat)
            mag = floor( log10( abs( maxprm(ii) ) ) );
            if mag <= 1
                ax.YAxis.Exponent = 0;
            else
                ax.YAxis.Exponent = mag;
            end


            set(gca,'LineWidth',2)
        else
            %   row < column -> upper triangle
            y = samps(:,ii);
            x = samps(:,jj);
            z = exp(samps(:,end));

            xbuffer = buffer*max( abs( [ minprm(jj),maxprm(jj) ] ) );
            ybuffer = buffer*max( abs( [ minprm(ii),maxprm(ii) ] ) );
            
            scatter(x,y,[],z,...
                "filled",'MarkerEdgeColor','k')
            hold on

             %   95% for x-axis
            xplot1 = quants(1,jj)*ones(1,100);
            xplot2 = quants(2,jj)*ones(1,100);
            yplot1 = linspace(minprm(ii)-ybuffer,maxprm(ii)+ybuffer);
            plot(xplot1,yplot1,...
                '--k','LineWidth',2)
            plot(xplot2,yplot1,...
                '--k','LineWidth',2)
            %   95% for y-axis
            yplot2 = quants(1,ii)*ones(1,100);
            yplot3 = quants(2,ii)*ones(1,100);
            xplot3 = linspace(minprm(jj)-xbuffer, maxprm(jj)+xbuffer);
            plot(xplot3, yplot2, ...
                '--k','LineWidth',2)
            plot(xplot3, yplot3, ...
                '--k','LineWidth',2)

            %   Show location of best sample
            xbest = samps(idx,jj);
            ybest = samps(idx,ii);
            plot(xbest*ones(1,100),yplot1,...
                "Color",[0.8500 0.3250 0.0980],...
                "LineWidth",2)
            plot(xplot3, ybest*ones(1,100),...
                "Color",[0.8500 0.3250 0.0980],...
                'LineWidth',2)
            scatter(xbest,ybest,50,[0.8500 0.3250 0.0980],...
                'filled','o',...
                'MarkerEdgeColor','k')

            box on
            hold off
            tickvec = linspace(minprm(jj), maxprm(jj), n_ticks+2);
            xticks(tickvec(2:end-1))

            tickvec = linspace(minprm(ii), maxprm(ii), n_ticks+2);
            yticks(tickvec(2:end-1))

            xticklabels('')
            yticklabels('')

            xlim([minprm(jj)-xbuffer, maxprm(jj)+xbuffer])
            ylim([minprm(ii)-ybuffer, maxprm(ii)+ybuffer])


            set(gca,'LineWidth',2)
        end
    end
end
%   Plot aesthetics
cb = colorbar('Ticks',[0.1, 0.9],...
    'TickLabels',{'Poor fit\newline{to data}', 'Good fit\newline{to data}'});
cb.Layout.Tile = 'east';
end