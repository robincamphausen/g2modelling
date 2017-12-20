function [transmitted, reflected] = mirror(input,transmissionProbability)
%mirror function: gives transmitted and reflected photon pulsetrains for a given
%input state and a specified transmissionProbability
% Author: Robin Camphausen
% Last modified: 19/12/2017

%   For each timebin the number of photons transmitted is determined by
%   drawing a number k from the binomial distribution with n trials and
%   probability p (according to the standard definition of the binomial
%   distribution), where n is the number of photons in that timebin (i.e.
%   the ith element of input) and p is transmissionProbability.
% -------------------------------------------------------------------------

% Check if input is pure Fock state - i.e. that n is the same for every
% input element:
if sum(abs(input-mean(input)))==0
    n = mean(input);
    transmitted = binornd(n,transmissionProbability,1,length(input));

% if n is not the same for every input element, calculate for every
% separate value of n:
else
%     fast method:
% tic
    transmitted = zeros(1, length(input));
    n = 0:max(input);
    sumOfSums = 0;
    for i = 1:length(n)
        calculateForThisN = input ==n(i);
        transmitted(calculateForThisN) = ...
            binornd(n(i),transmissionProbability,1,sum(calculateForThisN));
        sumOfSums = sumOfSums + sum(calculateForThisN);
    end
%     display(sumOfSums)

%     toc
%     tic
% %     original (stupid) method:
%     transmitted2 = zeros(1, length(input));
%     for i = 1:length(input)
%         transmitted2(i) = binornd(input(i),transmissionProbability,1);
%     end
%     toc
end
% display('Difference between original and new: ')
% sum(transmitted2) - sum(transmitted)
% reflected is just the difference between input and transmitted
reflected = input - transmitted;    
end

