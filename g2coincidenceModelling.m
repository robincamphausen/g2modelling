% g2 coincidence modelling script
% Author: Robin Camphausen
% Last modified: 19/12/2017

% Description: Takes one or several photon inputs (coherent or Fock states) and
% models photon loss (binomial dist. probability) and splitting through a
% mirror (also binomial dist.), giving 2-detector coincidence at the end of
% this whole process.
% -------------------------------------------------------------------------
decays = 4:0.5:9;
for loopyDecays = 1:length(decays)

clearvars -except decays loopyDecays
close all

% Specify input states:
numPulses = 1000000; %number of consecutive pulses passing through set-up per loop
pulsePeriod = 12.5; %in ns

% Number of loops - i.e. total number of photon states considered is
% numLoops*numPulses
numLoops = 100000;
% numLoops = 100;

N_Fock = 1; %number of photons per pulse if choosing Fock state
N_coherent = 10; %mean number of photons per pulse if choosing coherent state
whatTheFock1 = 1; %if ==1 then input1 is a Fock state, if ==0 no input1
whatTheFock2 = 0; %if ==1 then input2 is a Fock state, if ==0 no input2

% Specify input state decay lifetimes (exponential decay):
% tauDecay1 = 6.5; %in ns
tauDecay1 = decays(loopyDecays); %in ns
tauDecay2 = 5; %in ns

% Specify system loss and detector dead time:
transmissionProb = 0.001; %transmissionProb==1 means no loss, ==0 means perfect loss
deadTime = 2000; % in ns

% timebin for plotting:
bin = 1;

% Initialise coincidences and time variables:
coincidences = [];
timeInit = 0;
timeLoss = 0;
timeBS = 0;
timeDT = 0;
timeSortingCoinc = 0;

for loopLooper = 1:numLoops

% -------------------------------------------------------------------------
% Initialise pulsetrain variables and apply exponential decays according to
% the input state decay lifetimes
tic

% Inputs are 2*numPulses matrices, where the first row gives the number of
% photons in each time bin and the second row the timing of each bin:
if whatTheFock1 ~= 0
    input1 = zeros(2,numPulses);
    input1(2,:) = 0:pulsePeriod:(pulsePeriod*numPulses - pulsePeriod);
% Row 2 initialises each timebin separated by the pulsePeriod
    if whatTheFock1 == 1
        input1(1,:) = N_Fock;
% If Fock state, number of photons is the same for each timebin.
    else
        input1(1,:) = poissrnd(N_coherent, 1, numPulses);
% If coherent state, number of photons given by Poisson distribution with
% mean N_coherent
    end
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
end

% Apply decay lifetimes and merge inputs:
if whatTheFock1 ==0 && whatTheFock2 ==0
    display('No inputs!')
    return
    
% Add time randomly drawn from exponential decay with decay lifetime
% tauDecay1 or tauDecay2 respectively:
elseif whatTheFock1 ~=0 && whatTheFock2 ==0
    input1(2,:) = input1(2,:) - tauDecay1*log(rand(1,numPulses));
    input = input1;
elseif whatTheFock1 ==0 && whatTheFock2 ~=0
    input2(2,:) = input2(2,:) - tauDecay2*log(rand(1,numPulses));
    input = input2;
else
% If we have both input1 and input2, sort accordding to timestamps:
    input1(2,:) = input1(2,:) - tauDecay1*log(rand(1,numPulses));
    input2(2,:) = input2(2,:) - tauDecay2*log(rand(1,numPulses));
    input = [input1 input2];
    input = sortrows(input',2)';
end
clear input1 input2
timeInit = timeInit + toc;

% -------------------------------------------------------------------------
tic
% Pass through loss element and through mirror, and apply detections of photons:

afterLoss = mirror(input(1,:),transmissionProb);
timeLoss = timeLoss + toc;
tic
[mirrorT, mirrorR] = mirror(afterLoss,0.5); %mirror has 50% transmission prob
mirrorT(2,:) = input(2,:); %transmitted photons get their timestamps
mirrorR(2,:) = input(2,:); %reflected photons get their timestamps
timeBS = timeBS + toc;
% Apply detection: Call Transmitted detections channel 0 and Reflected
% detections channel 1
tic
ch0_detected = mirrorT(1,:) ~=0;
ch0_timestamps = mirrorT(2,ch0_detected);
ch0_timestamps = deadtime(ch0_timestamps,deadTime); %take into account DT
ch0_timestamps(2,:) = -0.5; %channel marker (as ch0-ch1 is positive delta tau)

ch1_detected = mirrorR(1,:) ~=0;
ch1_timestamps = mirrorR(2,ch1_detected);
ch1_timestamps = deadtime(ch1_timestamps,deadTime); %take into account DT
ch1_timestamps(2,:) = 0.5; %channel marker (as ch1-ch0 is negative delta tau)

timeDT = timeDT + toc;
% -------------------------------------------------------------------------
tic
% Calculate coincidences:
% If ith ch0 was detected before ch1 (i.e. ch0_timestamp < ch1_timestamp) the
% channel marker will be 0.5-(-0.5) = 1, thus interphotonDelays(2,i) =1.
% Likewise if ch1 was detected before ch0, will have -0.5-0.5 = -1.

timestamps = [ch0_timestamps ch1_timestamps]; %merge timestamp arrays
timestamps = sortrows(timestamps',1)'; %order total timestamp array in ascending order
interphotonDelays = diff(timestamps,1,2);
% Disregard interphoton delays between consecutive same channel detections:
sameChannelDelays = interphotonDelays(2,:)==0;
interphotonDelays(:,sameChannelDelays) = []; %this removes same channel coincidences

coincidencesFromThisLoop = interphotonDelays(1,:).*interphotonDelays(2,:);
% display('sorting and coincidences:')
timeSortingCoinc = timeSortingCoinc + toc;

coincidences = [coincidences, coincidencesFromThisLoop];

clear afterLoss ch0_detected ch1_detected ch0_timestamps ch1_timestamps ...
    input mirrorR mirrorT
end

% -------------------------------------------------------------------------
tic
% Plot histogram of coincidences:
figure;
lbx = 'Relative Delay [ns]';
lby = 'Coincidences';
h1 = histogram(coincidences, -100: bin :100);

h1.FaceColor = [0.5, 0.5, 0.5];
fontsize = 20;

set(gca,'FontSize',24);
% set(gcf, 'Position', [100, 300, 900, 700]);
pbaspect([9 7 1]);
display('plotting:')
toc
display(timeInit)
display(timeLoss)
display(timeBS)
display(timeDT)
display(timeSortingCoinc)

filename = strcat('singlePhoton_tauDecay1_', num2str(tauDecay1), 'ns.mat');
save(filename)
end