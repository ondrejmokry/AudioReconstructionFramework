function [ SDRs, SDRs_deq, SDRs_inp, ODGs, times ] = meaner(fold,gflag)
% MEANER Compute the average performance of the reconstruction procedure
% based on several signals.
%
% Input arguments
%       fold ... path to the folder with the results, ending with '/'
%       gflag .. logical switch between the two variants of the assignment,
%                for the results using the function g (condatg.m,
%                multitestg.m), set true, otherwise (condat.m, multitest.m)
%                set false
%
% Output arguments
%       SDRs, SDRs_deq, SDRs_inp, ODGs, times ..... mean of the performance
%                 measures and execution times; for the structure of the
%                 arrays, see the functions multitest.m or multitestg.m
%
% Date: 21/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% initialization
ids      = {'a08','a16','a18','a25','a35','a41','a42','a58','a60','a66'};

% load the first file to read sizes of the arrays
L        = load([ fold, 'experiment_a08.mat' ]);
SDRs     = zeros(size(L.SDRs));
SDRs_deq = zeros(size(L.SDRs_deq));
SDRs_inp = zeros(size(L.SDRs_inp));
ODGs     = zeros(size(L.ODGs));
times    = zeros(size(L.times));

%% cycle through the signals (and flames)
for n = 1:10
    if gflag
        filename = [ fold, 'experimentg_', ids{n}, '.mat' ];
    else
        filename = [ fold, 'experiment_', ids{n}, '.mat' ];
    end
    if ~isfile(filename)
        break
    else
        L = load(filename);
        SDRs     = SDRs     + L.SDRs;
        SDRs_deq = SDRs_deq + L.SDRs_deq;
        SDRs_inp = SDRs_inp + L.SDRs_inp;
        ODGs     = ODGs     + L.ODGs;
        times    = times    + L.times;
    end 
end

% correction of the break situation
if n < 10
    n = n - 1;
end

% divide by number of signals
SDRs     = SDRs/n;
SDRs_deq = SDRs_deq/n;
SDRs_inp = SDRs_inp/n;
ODGs     = ODGs/n;
times    = times/n;

end