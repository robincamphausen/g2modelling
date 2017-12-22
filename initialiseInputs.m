function [ input ] = ...
    initialiseInputs( whatTheFock1,whatTheFock2,numPulses,pulsePeriod,...
    N_Fock, N_coherent, tauDecay1, tauDecay2 )
% Initialise pulsetrain variables and apply exponential decays according to
% the input state decay lifetimes
% Author:Robin Camphausen
% Last Modified: 22/12/2017
% ------------------------------------------------------------------------
input = [];

if whatTheFock1 ==0 && whatTheFock2 ==0
    display('No inputs!')
    return
end

% Inputs are 2*numPulses matrices, where the first row gives the number of
% photons in each time bin and the second row the timing of each bin:
if whatTheFock1 ~= 0
    input = zeros(2,numPulses);
    input(2,:) = 0:pulsePeriod:(pulsePeriod*numPulses - pulsePeriod);
% Row 2 initialises each timebin separated by the pulsePeriod
    if whatTheFock1 == 1
        input(1,:) = N_Fock;
% If Fock state, number of photons is the same for each timebin.
    else
        input(1,:) = poissrnd(N_coherent, 1, numPulses);
% If coherent state, number of photons given by Poisson distribution with
% mean N_coherent
    end

% get rid of zero elements (if coherent state has timebins with no photons)
inputHasNoPhotons = input(1,:)==0;
input(:,inputHasNoPhotons) = [];
% Add time randomly drawn from exponential decay with decay lifetime
% tauDecay1:
input(2,:) = input(2,:) - tauDecay1*log(rand(1,length(input(1,:))));
end


if whatTheFock2 ~= 0
    input2 = zeros(2,numPulses);
    input2(2,:) = 0:pulsePeriod:(pulsePeriod*numPulses - pulsePeriod);
% Row 2 initialises each timebin separated by the pulsePeriod
    if whatTheFock2 == 1
        input2(1,:) = N_Fock;
% If Fock state, number of photons is the same for each timebin.
    else
        input2(1,:) = poissrnd(N_coherent, 1, numPulses);
% If coherent state, number of photons given by Poisson distribution with
% mean N_coherent
    end

% get rid of zero elements (if coherent state has timebins with no photons)
input2HasNoPhotons = input2(1,:)==0;
input2(:,input2HasNoPhotons) = [];
% Add time randomly drawn from exponential decay with decay lifetime
% tauDecay2:
input2(2,:) = input2(2,:) - tauDecay2*log(rand(1,length(input2(1,:))));

input = [input input2];
end
% maybe don't need to sort here, as we are sorting later anyway?
% input = sortrows(input',2)';

% if whatTheFock1 ==0 && whatTheFock2 ==0
%     display('No inputs!')
%     return
%     
% % Add time randomly drawn from exponential decay with decay lifetime
% % tauDecay1 or tauDecay2 respectively:
% elseif whatTheFock1 ~=0 && whatTheFock2 ==0
%     input1(2,:) = input1(2,:) - tauDecay1*log(rand(1,numPulses));
%     input = input1;
% elseif whatTheFock1 ==0 && whatTheFock2 ~=0
%     input2(2,:) = input2(2,:) - tauDecay2*log(rand(1,numPulses));
%     input = input2;
% else
% 
%     
% % If we have both input1 and input2, sort accordding to timestamps:
%     input1(2,:) = input1(2,:) - tauDecay1*log(rand(1,length(input1(1,:))));
%     input2(2,:) = input2(2,:) - tauDecay2*log(rand(1,length(input2(1,:))));
%     input = [input1 input2];
%     input = sortrows(input',2)';
% end

end

