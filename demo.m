%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           demo file for the general reconstruction framework            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date: 29/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

clear
clc
close all

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

fprintf('\nAvailable signals:\n')
for i = 1:length(sigs)
    fprintf('  %2d.  %s\n',i,sigs{i})
end

prompt = '\nChoose signal number (1-10): ';
signum = input(prompt);

prompt = '\nChoose bit depth (1-32): ';
nbits  = input(prompt);

prompt = '\nChoose percentage of reliable samples (0-100): ';
pT     = input(prompt)/100;

prompt = '\nChoose percentage of reliable coefficients (0-100): ';
pTF    = input(prompt)/100;

prompt = '\nSave wavs? (0/1): ';
wavs   = input(prompt);

fprintf('\nPlease wait for a while...')

%% take a signal
signame = sigs{signum};
signal  = eval(signame);
signal  = signal/max(abs(signal));

%% choose TF transform
F = frametight(frame('dgt',{'sine',2048},1024,2048,'timeinv'));

%% resynthesize the signal so that it has the right length
signal = real(frsyn(F, frana(F,signal)));

%% precompute the random vector for the sake of reproducibility
rands = rand(size(signal));

%% parameters
wT  = nbits; % wordlength in T domain
wTF = nbits; % wordlength in TF domain

%% quantize the coefficients and the signal
c = frana(F, signal);
[cq, dTF] = cquant(c, wTF);
[sq, dT]  = quant(signal, wT);

%% generate the T domain mask
maskT = rands <= pT;

%% generate the TF domain mask
% we take the coefficients that are largest in magnitude
crel   = hard_thresholding(c, floor(pTF*length(c)/2));
maskTF = logical(abs(crel));

%% drop coefficients and signal samples
sq = sq.*maskT;
cq = cq.*maskTF;

%% set the parameters of the Condat algorithm
% parameters of the model
model.ana    = @(x) frana(F, x);
model.syn    = @(x) frsyn(F, x);
model.projT  = @(x) x.*(~maskT) + proj(x, sq, dT).*maskT;
model.projTF = @(x) x.*(~maskTF) + cproj(x, cq, dTF).*maskTF;
model.dim    = [ length(signal), length(c) ];

% parameters of the algorithm
algo.tau   = sqrt(1/2);
algo.sigma = sqrt(1/2);
algo.rho   = 1;
algo.maxit = 300;
algo.tol   = 0;

% sparsifying step
model.sparse = @(x) sign(x) .* max(abs(x) - 1/algo.sigma, 0);

%% run the Condat algorithm
tic
[xana, relnormana] = condatg('analysis', model, algo);
xana = real(xana);
t = toc;
fprintf('\nAnalysis model done, time %0.2f s.\n',t)
fprintf('   SDR on missing samples:   %0.2f dB\n', sdr(signal(~maskT), xana(~maskT)));
fprintf('   SDR on quantized samples: %0.2f dB\n', sdr(signal(maskT),  xana(maskT)));
fprintf('   SDR on the whole signal:  %0.2f dB\n', sdr(signal, xana));

tic
[csyn, relnormsyn] = condatg('synthesis', model, algo);
xsyn = real(frsyn(F, csyn));
t = toc;
fprintf('\nSynthesis model done, time %0.2f s.\n',t)
fprintf('   SDR on missing samples:   %0.2f dB\n', sdr(signal(~maskT), xsyn(~maskT)));
fprintf('   SDR on quantized samples: %0.2f dB\n', sdr(signal(maskT),  xsyn(maskT)));
fprintf('   SDR on the whole signal:  %0.2f dB\n', sdr(signal, xsyn));

figure
hold on
semilogy(relnormana)
semilogy(relnormsyn)
legend('analysis model','synthesis model')
xlabel('iteration number')
ylabel('relative norm of the solution')
title('Convergence plot')

%% plot the observed signal and coefficients
limits = [2 2.01];

sq(~maskT) = NaN;

figure
subplot(1,3,1)
hold on
stem((1:length(signal))/fs,signal,'rx','linestyle','none','markersize',6);
stem((1:length(signal))/fs,sq,'bo','linestyle','none','markersize',6);
legend('original','observed','location','southeast')
title('T domain, observed (detail)')
grid on
set(gca,'XMinorGrid','on')
set(gca,'YMinorGrid','on')
box on
xlim(limits)
xlabel('Time (s)')

subplot(1,3,2)
plotframe(F,c,'dynrange',80,'fs',fs)
title('TF domain, original')
ylim([0 inf])

subplot(1,3,3)
plotframe(F,cq,'dynrange',80,'fs',fs)
title('TF domain, observed')
ylim([0 inf])

%% plot the results
figure
subplot(1,3,1)
hold on
plot((1:length(signal))/fs,signal)
plot((1:length(signal))/fs,xana)
plot((1:length(signal))/fs,xsyn)
legend('original','analysis model','synthesis model')
title('T domain, restored')
grid on
set(gca,'XMinorGrid','on')
set(gca,'YMinorGrid','on')
box on
xlabel('Time (s)')

subplot(1,3,2)
plotframe(F,frana(F,xana),'dynrange',80,'fs',fs)
title('TF domain, restored, analysis')
ylim([0 inf])

subplot(1,3,3)
plotframe(F,csyn,'dynrange',80,'fs',fs)
title('TF domain, restored, synthesis')
ylim([0 inf])

%% wavs
if wavs
    audiowrite(['results/wavs/',sigs{signum}, '_reference.wav'],signal,fs);
    audiowrite(['results/wavs/',sigs{signum}, '_analysis.wav'] ,xana,fs);
    audiowrite(['results/wavs/',sigs{signum}, '_synthesis.wav'],xsyn,fs);
end