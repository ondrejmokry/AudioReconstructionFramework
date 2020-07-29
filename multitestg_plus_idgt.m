%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           test of the degradation (amd reconstruction) model            %
%       with both T and TF domains quantized and partially missing        %
%            (using the function g of the Condat-Vu algorithm)            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date: 21/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

ltfatstart
rng(0)
addpath(genpath('PemoQ'));

% load signals
load('signals/EBU_SQAM.mat');
sigs = { 'a08_violin',...
         'a16_clarinet',...
         'a18_bassoon',...
         'a25_harp',...
         'a35_glockenspiel',...
         'a41_celesta',...
         'a42_accordion',...
         'a58_guitar_sarasate',...
         'a60_piano_schubert',...
         'a66_wind_ensemble_stravinsky' };

% set parameters
nbits = [ 2 4 8 16 32 ];
pTs   = 0.1 : 0.1 : 0.9;
pTFs  = 0.1 : 0.1 : 0.9;

% timer
t_start = clock;

combinationcounter = 0;
for signum = 1:length(sigs)
    
    ODGs      = NaN(7, length(nbits), length(pTs), length(pTFs));
    SDRs_inp  = NaN(7, length(nbits), length(pTs), length(pTFs));
    SDRs_deq  = NaN(7, length(nbits), length(pTs), length(pTFs));
    SDRs      = NaN(7, length(nbits), length(pTs), length(pTFs));
    times     = NaN(7, length(nbits), length(pTs), length(pTFs));

    %% take a signal and make it shorter
    signame = sigs{signum};
    s       = eval(signame);
    s       = s(fs+1:2*fs);
    s       = s/max(abs(s));

    L = load(['results/experiments/experimentg_', signame(1:3), '.mat']);
    ODGs    (1:6,:,:,:) = L.ODGs;
    SDRs_inp(1:6,:,:,:) = L.SDRs_inp;
    SDRs_deq(1:6,:,:,:) = L.SDRs_deq;
    SDRs    (1:6,:,:,:) = L.SDRs;
    times   (1:6,:,:,:) = L.times;

    %% curl the ends
    cosine   = cos(linspace(-pi/2,pi/2,400)').^2;
    s(1:200) = s(1:200).*cosine(1:200);
    s(end-199:end) = s(end-199:end).*cosine(201:end);

    %% choose TF transform
    F = frametight(frame('dgt',{'sine',2048},1024,2048,'timeinv'));

    %% resynthesize the signal so that it has the right length
    s = real(frsyn(F, frana(F,s)));

    %% precompute the random vector for the sake of reproducibility
    rands = rand(size(s));

    for i = 1:length(nbits)
        for j = 1:length(pTs)
            for k = 1:length(pTFs)

                combinationcounter = combinationcounter + 1;

                %% parameters
                wT  = nbits(i); % wordlength in T domain
                pT  = pTs(j);   % percentage of T domain samples available
                wTF = nbits(i); % wordlength in TF domain
                pTF = pTFs(k);  % percentage of TF coefficients available

                %% quantize the coefficients and the signal
                c = frana(F, s);
                [cq, dTF] = cquant(c, wTF);
                [sq, dT]  = quant(s, wT);

                %% generate the T domain mask
                maskT = rands <= pT;

                %% generate the TF domain mask
                % we take the coefficients that are largest in magnitude
                crel = hard_thresholding(c, floor(pTF*length(c)/2));
                maskTF = logical(abs(crel));

                %% drop coefficients and signal samples
                sq = sq.*maskT;
                cq = cq.*maskTF;

                %% set the parameters of the Condat algorithm
                model.ana = @(x) frana(F, x);
                model.syn = @(x) frsyn(F, x);
                model.projT  = @(x) x.*(~maskT) + proj(x,  sq, dT).*maskT;
                model.projTF = @(x) x.*(~maskTF) + cproj(x, cq, dTF).*maskTF;
                model.dim = [ length(s), length(c) ];

                algo.tau   = sqrt(1/2);
                algo.sigma = sqrt(1/2);
                algo.rho   = 1;
                algo.maxit = 300;
                algo.tol   = 0;

                model.sparse = @(x) sign(x) .* max(abs(x) - 1/algo.sigma, 0);
                
                %% compute simple idgt of the observed coefficients
                tic
                xidgt = real(frsyn(F,cq));
                times(7, i, j, k) = toc;

                %% compute the metrics
                fprintf('\nnbits: %d, pT: %.1f, pTF: %.1f\n', wT, pT, pTF)
                [~, ~, ODG, ~]       = audioqual(s, xidgt, fs);
                ODGs(7, i, j, k)     = ODG;
                SDRs_inp(7, i, j, k) = sdr(s(~maskT), xidgt(~maskT));
                SDRs_deq(7, i, j, k) = sdr(s(maskT), xidgt(maskT));
                SDRs(7, i, j, k)     = sdr(s, xidgt);

                save(['results/experiments_plus_idgt/experimentg_', signame(1:3), '.mat'],'ODGs','SDRs_inp','SDRs_deq','SDRs','times','nbits','pTs','pTFs','s','fs')
                
                % timer again
                t_now = clock;
                fprintf('\nSo far, the experiment has taken %d hours.',round(etime(t_now,t_start)/3600))
                estimatedtotalhours =...
                    etime(t_now,t_start)*length(sigs)*length(nbits)*length(pTs)*length(pTFs)...
                    /( combinationcounter * 3600);
                fprintf('\nEstimated remaining time: %d hours.\n',round(estimatedtotalhours - etime(t_now,t_start)/3600))
            end
        end
    end
end