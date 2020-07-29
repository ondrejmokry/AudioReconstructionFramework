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
    
    ODGs      = NaN(6, length(nbits), length(pTs), length(pTFs));
    SDRs_inp  = NaN(6, length(nbits), length(pTs), length(pTFs));
    SDRs_deq  = NaN(6, length(nbits), length(pTs), length(pTFs));
    SDRs      = NaN(6, length(nbits), length(pTs), length(pTFs));
    times     = NaN(6, length(nbits), length(pTs), length(pTFs));

    %% take a signal and make it shorter
    signame = sigs{signum};
    s       = eval(signame);
    s       = s(fs+1:2*fs);
    s       = s/max(abs(s));

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

                %% run the Condat algorithm
                tic
                [xana, ~] = condatg('analysis', model, algo);
                xana      = real(xana);
                times(1, i, j, k) = toc;

                tic
                [csyn, ~] = condatg('synthesis', model, algo);
                xsyn      = real(frsyn(F, csyn));
                times(2, i, j, k) = toc;

                %% run inpainting + dequentization in T domain
                model.projT  = @(x) x.*(~maskT) + proj(x, sq, dT).*maskT;
                model.projTF = @(x) x;

                tic
                [xanainp, ~] = condatg('analysis', model, algo);
                xanainp      = real(xanainp);
                times(3, i, j, k) = toc;

                tic
                [csyninp, ~] = condatg('synthesis', model, algo);
                xsyninp      = real(frsyn(F, csyninp));
                times(4, i, j, k) = toc;

                %% run inpainting + dequentization in TF domain
                model.projT  = @(x) x;
                model.projTF = @(x) x.*(~maskTF) + cproj(x, cq, dTF).*maskTF;
                
                tic
                [xanafre, ~] = condatg('analysis', model, algo);
                xanafre = real(xanafre);
                times(5, i, j, k) = toc;
                
                tic
                [csynfre, ~] = condatg('synthesis', model, algo);
                xsynfre = real(frsyn(F, csynfre));
                times(6, i, j, k) = toc;

                %% compute the metrics
                fprintf('\nnbits: %d, pT: %.1f, pTF: %.1f\n', wT, pT, pTF)
                for m = 1:6
                    switch m
                        case 1
                            rec = xana;
                            fprintf('\nanalysis model, both domains:\n')
                        case 2
                            rec = xsyn;
                            fprintf('\nsynthesis model, both domains:\n')
                        case 3
                            rec = xanainp;
                            fprintf('\nanalysis model, time domain:\n')
                        case 4
                            rec = xsyninp;
                            fprintf('\nsynthesis model, time domain:\n')
                        case 5
                            rec = xanafre;
                            fprintf('\nanalysis model, TF domain:\n')
                        case 6
                            rec = xsynfre;
                            fprintf('\nsynthesis model, TF domain:\n')
                    end
                    [~, ~, ODG, ~]       = audioqual(s, rec, fs);
                    ODGs(m, i, j, k)     = ODG;
                    SDRs_inp(m, i, j, k) = sdr(s(~maskT), rec(~maskT));
                    SDRs_deq(m, i, j, k) = sdr(s(maskT), rec(maskT));
                    SDRs(m, i, j, k)     = sdr(s, rec);

                    %% output to command window
                    fprintf(repmat('\b', 1, 22))
                    fprintf('   ODG: %f\n', ODG)
                    fprintf('   SDR on missing samples: %f\n', sdr(s(~maskT), rec(~maskT)));
                    fprintf('   SDR on quantized samples: %f\n', sdr(s(maskT), rec(maskT)));
                    fprintf('   SDR on the whole signal: %f\n', sdr(s, rec));

                end           
                save(['results/experiments/experimentg_', signame(1:3), '.mat'],'ODGs','SDRs_inp','SDRs_deq','SDRs','times','nbits','pTs','pTFs','s','fs')
                
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