%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%         plotting file for the general reconstruction framework          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date: 29/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

close all

%% load data
nbits = [2 4 8 16 32];
pTs   = 0.1:0.1:0.9;
pTFs  = 0.1:0.1:0.9;
fs    = 44100;
[ SDRs, SDRs_deq, SDRs_inp, ODGs, times ] = meaner('experiments_plus_idgt/',true);
map   = 'parula';

plotbars  = true;
plotlines = true;
metrics   = [1 4];
% 1 ODG
% 2 SDR on missing samples
% 3 SDR on quantized samples
% 4 SDR on the whole signal

%% bar plots

% Compare the models for given set of bits.
% The TF representation is redundant with redundancy 2.
%
% The number of bits used by the double-domain model is
%
%       pT * nbits + 2 * pTF * nbits = (pT + 2 * pTF) * nbits,
%
% multiplied with the length of the signal, which is ommited now.
%
% The number of bits used by the T-domain-only model is
%
%       pT * nbits.
%
% For the comparison, we fix the number of bits used in the T-domain-only
% model and search for the possible settings in the double-domain model.
%
%  only T       both T and TF
%       
%  pT = 0.1     cannot be reached from this experiment
%  pT = 0.2     cannot be reached from this experiment
%  pT = 0.3     pT = 0.1, pTF = 0.1
%  pT = 0.4     pT = 0.2, pTF = 0.1
%  pT = 0.5     pT = 0.3, pTF = 0.1 or pT = 0.1, pTF = 0.2 
%  pT = 0.6     pT = 0.4, pTF = 0.1 or pT = 0.2, pTF = 0.2
%  pT = 0.7     pT = 0.5, pTF = 0.1 or pT = 0.3, pTF = 0.2 or pT = 0.1, pTF = 0.3
%  pT = 0.8     pT = 0.6, pTF = 0.1 or pT = 0.4, pTF = 0.2 or pT = 0.2, pTF = 0.3
%  pT = 0.9     pT = 0.7, pTF = 0.1 or pT = 0.5, pTF = 0.2 or pT = 0.3, pTF = 0.3 or pT = 0.1, pTF = 0.4

if plotbars
    colors = eval([map, '(5)']);
    for bitnum = 1:5
        figure
        colormap(map)
        for subnum = 1:4
            subplot(2,2,subnum)
            switch subnum
                case {1,3}
                    data = ODGs;
                    plottitle = 'ODG';
                case {2,4}
                    data = SDRs;
                    plottitle = 'SDR on the whole signal';
            end        
            ana = ceil(subnum/2);
            plotdata        = NaN(9,4);
            plotdata(1:9,1) = data(2+ana,bitnum,:,1);
            plotdata(3:9,2) = data(ana,bitnum,1:7,1);    
            plotdata(5:9,3) = data(ana,bitnum,1:5,2);    
            plotdata(7:9,4) = data(ana,bitnum,1:3,3);   
            plotdata(9,5)   = data(ana,bitnum,1,4);
            b = bar(plotdata);
            legend(...
                '$$p_\mathrm{TF} = 0\, \% $$',...
                '$$p_\mathrm{TF} = 10\,\% $$',...
                '$$p_\mathrm{TF} = 20\,\% $$',...
                '$$p_\mathrm{TF} = 30\,\% $$',...
                '$$p_\mathrm{TF} = 40\,\% $$',...
                'location','northwest')
            xticklabels({'10\,\%','20\,\%','30\,\%','40\,\%','50\,\%','60\,\%','70\,\%','80\,\%','90\,\%'})
            xlabel('pct.\ of the reliable information (wrt the T domain)')
            switch subnum
                case {1,3}
                    ylabel('ODG')
                case {2,4}
                    ylabel('SDR [dB]')
            end 
            if ana == 1
                title([plottitle, ', analysis model'])
            else
                title([plottitle, ', synthesis model'])
            end
            for i = 1:5
                b(i).FaceColor = colors(i,:);
            end
        end
        sgtitle(['nbits = ', num2str(nbits(bitnum))])
    end
end

%% scatter plots
for fignum = metrics
    figure
    colormap(map)
    switch fignum
        case 1
            data = ODGs;
            plottitle = 'ODG';
        case 2
            data = SDRs_inp;
            plottitle = 'SDR on missing samples';
        case 3
            data = SDRs_deq;
            plottitle = 'SDR on quantized samples';
        case 4
            data = SDRs;
            plottitle = 'SDR on the whole signal';
    end
    
    for subnum = 1:2
        
        % re-organize data for single domain
        % 1     2  3   4     5     6
        % nbits pT pTF value bitsT bitsTF
        X = NaN(length(nbits)*length(pTs),6);
        counter = 1;
        for ii = 1:length(nbits)
            for jj = 1:length(pTs)
                X(counter,1) = nbits(ii);
                X(counter,2) = pTs(jj);
                X(counter,3) = 0;
                X(counter,4) = data(2+subnum, ii, jj, 1);
                counter      = counter + 1;
            end
        end
        X(:,5) = round(X(:,1).*X(:,2),1);
        X(:,6) = 0;
        
        % re-organize data for double domain
        % 1     2  3   4     5     6
        % nbits pT pTF value bitsT bitsTF
        XX = NaN(length(nbits)*length(pTs)*length(pTFs),6);
        counter = 1;
        for ii = 1:length(nbits)
            for jj = 1:length(pTs)
                for kk = 1:length(pTFs)
                    XX(counter,1) = nbits(ii);
                    XX(counter,2) = pTs(jj);
                    XX(counter,3) = pTFs(kk);
                    XX(counter,4) = data(subnum, ii, jj, kk);
                    counter       = counter + 1;
                end
            end
        end
        XX(:,5) = round(XX(:,1).*XX(:,2),1);
        XX(:,6) = round(2*XX(:,1).*XX(:,3),1);
        
        % find unique numbers of bits in single domain
        bitsS = unique(X(:,5));
       
        % for each number of bits used, take only the maximum value
        newX = NaN(length(bitsS),6);
        for bb = 1:length(bitsS)
            rows       = X(:,5) == bitsS(bb);
            maxval     = max(X(rows,4));
            row        = X(:,4) == maxval;
            newX(bb,:) = X(find(row,1),:); 
        end
        X = newX;
        
        % find unique numbers of bits in double domain (i.e. all the unique
        % combinations of the two domains)
        bitsS = unique(XX(:,5));
        bitsD = unique(XX(:,6));
        
        % for each number of bits used, take only the maximum value
        newXX   = NaN(length(bitsS)*length(bitsD),6);
        counter = 1;
        for bb = 1:length(bitsS)
            for cc = 1:length(bitsD)
                if bitsS(bb) + bitsD(cc) >= 32
                    continue
                end
                rowsT  = XX(:,5) == bitsS(bb);
                rowsTF = XX(:,6) == bitsD(cc);
                rows   = logical(rowsT.*rowsTF);
                if sum(rows) > 0
                    maxval           = max(XX(rows,4));
                    row              = XX(:,4) == maxval;
                    newXX(counter,:) = XX(find(row,1),:);
                    counter          = counter + 1;
                end
            end
        end
        XX = newXX(1:counter-1,:);
        
        % computing the bit-depths per second
        X(:,5) = X(:,5)*fs;
        XX(:,5:6) = XX(:,5:6)*fs;
        
        % initialize the plot
        subplot(1,2,subnum)
        hold on
             
        % read the colormap for equibitals
        cmap = colormap; % values of the colormap
        clim = [min([X(:,4); XX(:,4)]), max([X(:,4); XX(:,4)])]; % min and max values definining the colormap
        % clim(1) ~ 0 ~ the first color in cmap
        % clim(2) ~ 1 ~ the last color in cmap
        
        % plot the equibitals
        wdth = @(bits) log2(bits)^2/5;
        for index = 1:length(X)
            xval = X(index,5);
            cval = X(index,4);
            
            % find the color in the colormap
            colorindex = 1 + round(...
                (length(cmap)-1)*(cval-clim(1))/(clim(2)-clim(1)));
            
            % plot the line
            t = linspace(0, xval, 300);
            plot(t, xval - t, 'color', cmap(colorindex,:),'linewidth',wdth(X(index,1)))
        end
        
        % scatter the values
        sz = @(bits) 12*sqrt(bits);
        scatter(XX(:,5),XX(:,6),sz(XX(:,1)),XX(:,4),'filled','markeredgecolor',0.4*[1 1 1]);
        colorbar
        
        % make the awesome legend
        for n = 1:length(nbits)
            
            h(n) = plot([0 1], [0 0],...
                'linestyle','none',...
                'color',0.4*[1 1 1],...
                'marker','o',...
                'markerfacecolor',0.4*[1 1 1],...
                'markeredgecolor','none',...
                'markersize',sqrt(sz(nbits(n))));
            
            hh(n) = plot([0 1], [0 0],...
                'linewidth',wdth(nbits(n)),...
                'color',0.4*[1 1 1]);
        end
        legend([h hh],'2 bits','4 bits','8 bits','16 bits','32 bits',...
            '2 bits','4 bits','8 bits','16 bits','32 bits',...
            'location','southeast','numcolumns',2)
        
        % axes
        xlim([pTs(1)*fs,1.1*nbits(end)*pTs(end)*fs])
        ylim([pTs(1)*fs,1.1*nbits(end)*pTs(end)*fs])
        
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        axis square
        
        xlabel('bit rate in T domain [bits/s]')
        ylabel('bit rate in TF domain [bits/s]')
        
        % title
        if subnum == 1
            title('analysis model')
        else
            title('synthesis model')
        end
               
    end
    sgtitle(plottitle)
end

%% line plots
if plotlines
    for fignum = metrics
        figure
        colormap(map)
        switch fignum
            case 1
                data = ODGs;
                plottitle = 'ODG';
                ylab = 'ODG';
            case 2
                data = SDRs_inp;
                plottitle = 'SDR on missing samples';
                ylab = 'SDR [dB]';
            case 3
                data = SDRs_deq;
                plottitle = 'SDR on quantized samples';
                ylab = 'SDR [dB]';
            case 4
                data = SDRs;
                plottitle = 'SDR on the whole signal';
                ylab = 'SDR [dB]';
        end

        for subnum = 1:2

            % re-organize data for single domain
            % 1             2
            % total bitrate value
            X = NaN(length(nbits)*length(pTs),2);  % only T domain
            Y = NaN(length(nbits)*length(pTFs),2); % only TF domain
            Z = NaN(length(nbits)*length(pTFs),2); % IDGT of the quantized coefficients
            counter = 1;
            for ii = 1:length(nbits)
                for jj = 1:length(pTs) % suppose pTs = pTFs
                    X(counter,1) = nbits(ii)*pTs(jj);
                    Y(counter,1) = 2*nbits(ii)*pTs(jj);
                    Z(counter,1) = 2*nbits(ii)*pTs(jj);
                    X(counter,2) = data(2+subnum, ii, jj, 1);
                    Y(counter,2) = data(4+subnum, ii, 1, jj);
                    Z(counter,2) = data(7, ii, 1, jj);
                    counter      = counter + 1;
                end
            end
            X(:,1) = round(X(:,1),1);
            Y(:,1) = round(Y(:,1),1);
            Z(:,1) = round(Z(:,1),1);

            % re-organize data for double domain
            % 1             2
            % total bitrate value
            XX = NaN(length(nbits)*length(pTs)*length(pTFs),2);
            counter = 1;
            for ii = 1:length(nbits)
                for jj = 1:length(pTs)
                    for kk = 1:length(pTFs)
                        XX(counter,1) = nbits(ii)*pTs(jj) + 2*nbits(ii)*pTFs(kk);
                        XX(counter,2) = data(subnum, ii, jj, kk);
                        counter       = counter + 1;
                    end
                end
            end
            XX(:,1) = round(XX(:,1),1);

            % find unique numbers of bits in single domain
            bitsS = unique(X(:,1)); % S for single domain

            % for each number of bits used, take only the maximum value
            % up to the given number of bits
            newX = NaN(length(bitsS),2);
            for bb = 1:length(bitsS)
                rows       = X(:,1) <= bitsS(bb);
                newX(bb,1) = bitsS(bb);
                newX(bb,2) = max(X(rows,2));
            end
            X = newX;
            
            % do the same using only the TF domain and using the simple
            % IDGT of the quantized coefficients
            bitsS = unique(Y(:,1)); % S for single domain
            newY  = NaN(length(bitsS),2);
            newZ  = NaN(length(bitsS),2);
            counter = 1;
            for bb = 1:length(bitsS)
                if bitsS(bb) >= 32
                    continue
                end
                rows       = Y(:,1) <= bitsS(bb);
                newY(bb,1) = bitsS(bb);
                newY(bb,2) = max(Y(rows,2));
                newZ(bb,1) = bitsS(bb);
                newZ(bb,2) = max(Z(rows,2));
                counter    = counter + 1;
            end
            Y = newY(1:counter-1,:);
            Z = newZ(1:counter-1,:);

            % find unique numbers of bits in double domain (i.e. all the unique
            % combinations of the two domains)
            bitsD = unique(XX(:,1)); % D for double domain

            % for each number of bits used, take only the maximum value
            % up to the given number of bits
            newXX   = NaN(length(bitsD),2);
            counter = 1;
            for bb = 1:length(bitsD)
                if bitsD(bb) >= 32
                    continue
                end
                rows        = XX(:,1) <= bitsD(bb);
                newXX(bb,1) = bitsD(bb);
                newXX(bb,2) = max(XX(rows,2));
                counter     = counter + 1;
            end
            XX = newXX(1:counter-1,:);

            % computing the bit-depths per second
            X(:,1)  = X(:,1)*fs;
            Y(:,1)  = Y(:,1)*fs;
            Z(:,1)  = Z(:,1)*fs;
            XX(:,1) = XX(:,1)*fs;

            % initialize the plot
            subplot(1,2,subnum)
            hold on

            % plot the lines
            plot(X(:,1),X(:,2))
            plot(Y(:,1),Y(:,2))
            plot(Z(:,1),Z(:,2))
            plot(XX(:,1),XX(:,2))
            
            % legend
            legend('T-domain-only',...
                'TF-domain-only',...
                'direct synthesis',...
                'double-domain',...
                'location','northwest')

            % axes    
            set(gca,'xscale','log')
            xlabel('total bit rate [bits/s]')
            ylabel(ylab)
            xlim([min(XX(:,1)), Inf])
            if fignum > 1
                ylim([0 Inf])
            end
            grid on
            set(gca,'XMinorGrid','on')
            set(gca,'YMinorGrid','on')
            box on

            % title
            if subnum == 1
                title('analysis model')
            else
                title('synthesis model')
            end

        end
        sgtitle(plottitle)
    end
end