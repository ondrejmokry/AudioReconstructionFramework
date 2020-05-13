% (1) quantized coefficients in TF domain
% (2) some of them missing
% (3) quantized samples in T domain
% (4) some of them missing

PQ = exist('PemoQ','dir');
if PQ
    addpath(genpath('PemoQ'));
end

nbits = [ 2 4 8 16 32 ];
pTs   = 0.1 : 0.1 : 0.9;
pTFs  = 0.1 : 0.1 : 0.9;

for signal = 1:4
    
    switch signal
        case 1
            signame = 'violin';
        case 2
            signame = 'group_of_two';
        case 3
            signame = 'group_of_four';
        case 4
            signame = 'group_of_all';
    end

    ODGs      = NaN(4, length(nbits), length(pTs), length(pTFs));
    SDRs_inp  = NaN(4, length(nbits), length(pTs), length(pTFs));
    SDRs_deq  = NaN(4, length(nbits), length(pTs), length(pTFs));
    SDRs      = NaN(4, length(nbits), length(pTs), length(pTFs));

    %% take a signal
    [ s, fs ] = audioread(['signals/',signame,'.wav']);
    s = s(fs+1:2*fs);
    s = s/max(abs(s));

    %% choose TF transform
    F = frametight(frame('dgt',{'sine',2048},1024,2048,'timeinv'));

    %% resynthesize the signal so that it has the right length
    s = real(frsyn(F, frana(F,s)));

    %% precompute the random vector for the sake of reproducibility
    rands = rand(size(s));

    for i = 1:length(nbits)
        for j = 1:length(pTs)
            for k = 1:length(pTFs)

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

                algo.tau   = 0.5;
                algo.sigma = 0.5;
                algo.rho   = 1;
                algo.maxit = 300;
                algo.tol   = 0;
                
                model.sparse = @(x) sign(x) .* max(abs(x) - 1/algo.sigma, 0);

                %% run the Condat algorithm
                [xana, relnormana] = condat('analysis', model, algo);
                [csyn, relnormsyn] = condat('synthesis', model, algo);

                xana = real(xana);
                xsyn = real(frsyn(F, csyn));

                %% run inpainting + dequentization in T domain for reference
                [sqq, dTT] = quant(s, wT);
                model.projT  = @(x) x.*(~maskT) + proj(x, sqq, dTT).*maskT;
                model.projTF = @(x) x;

                [xanainp, ~] = condat('analysis', model, algo);
                [csyninp, ~] = condat('synthesis', model, algo);

                xanainp = real(xanainp);
                xsyninp = real(frsyn(F, csyninp));

                %% compute the metrics
                fprintf('\nnbits: %d, pT: %.1f, pTF: %.1f\n', wT, pT, pTF)           
                for l = 1:4
                    switch l
                        case 1
                            rec = xana;
                            fprintf('\nanalysis model, both domains:\n\n')
                        case 2
                            rec = xsyn;
                            fprintf('\nsynthesis model, both domains:\n\n')
                        case 3
                            rec = xanainp;
                            fprintf('\nanalysis model, time domain:\n\n')
                        case 4
                            rec = xsyninp;
                            fprintf('\nsynthesis model, time domain:\n\n')
                    end
                    if PQ
                        [~, ~, ODG, ~]   = audioqual(s/10, rec/10, fs);
                        ODGs(l, i, j, k) = ODG;
                    end
                    SDRs_inp(l, i, j, k) = sdr(s(~maskT), rec(~maskT));
                    SDRs_deq(l, i, j, k) = sdr(s(maskT), rec(maskT));
                    SDRs(l, i, j, k)     = sdr(s, rec);

                    %% output to command window
                    fprintf(repmat('\b', 1, 22))                
                    if PQ
                        fprintf('ODG: %f\n', ODG)
                    end
                    fprintf('SDR on missing samples: %f\n', sdr(s(~maskT), rec(~maskT)));
                    fprintf('SDR on quantized samples: %f\n', sdr(s(maskT), rec(maskT)));
                    fprintf('SDR on the whole signal: %f\n', sdr(s, rec));

                end           
                % save(['experiment_', signame, '.mat'],'ODGs','SDRs_inp','SDRs_deq','SDRs','nbits','pTs','pTFs','s','fs')
            end
        end
    end
end