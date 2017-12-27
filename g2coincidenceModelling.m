% g2 coincidence modelling script
% Author: Robin Camphausen
% Last modified: 22/12/2017

% Description: Takes one or several photon inputs (coherent or Fock states) and
% models photon loss (binomial dist. probability) and splitting through a
% mirror (also binomial dist.), giving 2-detector coincidence at the end of
% this whole process.
% -------------------------------------------------------------------------
decays1 = 1:3;
decays2 = 18:2:24;
for loopyDecays1 = 1:length(decays1)
    for loopyDecays2  = 1:length(decays2)

clearvars -except decays1 loopyDecays1 decays2 loopyDecays2
% clear all
close all

% Specify input states:
numPulses = 1e7; %number of consecutive pulses passing through set-up per loop
pulsePeriod = 12.5; %in ns

% Number of loops - i.e. total number of photon states considered is
% numLoops*numPulses
% numLoops = 1;
numLoops = 2000;

N_Fock = 1; %number of photons per pulse if choosing Fock state
N_coherent = 1; %mean number of photons per pulse if choosing coherent state

% if whatTheFock ==1 then input is a Fock state, if ==0 no input, if >=3 is coherent
% (make sure it's 3 or larger, due to the loopLooper condition, as
% whatTheFock1+whatTheFock2 < 3 if neither state is coherent - and thus
% there is no need to rerun initialisation of pulsetrain):
whatTheFock1 = 1; 
whatTheFock2 = 1;

% Specify input state decay lifetimes (exponential decay):
% tauDecay1 = 6.5; %in ns
tauDecay1 = decays1(loopyDecays1); %in ns
tauDecay2 = 5; %in ns
tauDecay2 = decays2(loopyDecays2); %in ns

% Specify system loss and detector dead time:
transmissionProb = 0.001; %transmissionProb==1 means no loss, ==0 means perfect loss
deadTime = 2000; % in ns

% timebin for plotting:
bin = 0.5;
plottingRange = 100;

% -------------------------------------------------------------------------
% Initialise pulsetrain variables and apply exponential decays according to
% the input state decay lifetimes

% Initialise coincidences and time variables:
coincidences = [];
timeInit = 0;
timeLoss = 0;
timeBS = 0;
timeDT = 0;
timeSortingCoinc = 0;
tic

% Photon input is 2*numPulses matrices, where the first row gives the number of
% photons in each time bin and the second row the timing of each bin (note
% that timebins with zero photons are always removed to speed up calculation)
input  = initialiseInputs( whatTheFock1,whatTheFock2,numPulses,...
    pulsePeriod, N_Fock, N_coherent, tauDecay1, tauDecay2);
% -------------------------------------------------------------------------

for loopLooper = 1:numLoops
    if loopLooper~=1 && (whatTheFock1+whatTheFock2)>=3
        input  = initialiseInputs( whatTheFock1,whatTheFock2,numPulses,...
    pulsePeriod, N_Fock, N_coherent, tauDecay1, tauDecay2);
% i.e. if one or both inputs is coherent, have to redo initialisation as
% coherent distribution is probabilistic, whereas Fock state is always the
% same
    end
timeInit = timeInit + toc;
% -------------------------------------------------------------------------
tic
% Pass through loss element and through mirror, and apply detections of photons:
pulsetrain(1,:) = mirror(input(1,:),transmissionProb);
pulsetrain(2,:) = input(2,:);
% remove timebins with zero photons in them again:
thisTimebinHasNoPhotons = pulsetrain(1,:)==0;
pulsetrain(:,thisTimebinHasNoPhotons) = [];
timeLoss = timeLoss + toc;

tic
% Beamsplitter and detection - Call Transmitted detections channel 0 and Reflected
% detections channel 1:
[ch0, ch1] = mirror(pulsetrain(1,:),0.5); %mirror has 50% transmission prob
ch0(2,:) = pulsetrain(2,:); %transmitted photons get their timestamps
ch1(2,:) = pulsetrain(2,:); %reflected photons get their timestamps

% once again remove timebins with no photons in them
thisTimebinHasPhotons = ch0(1,:)~=0;
ch0 = ch0(2,thisTimebinHasPhotons);
thisTimebinHasPhotons = ch1(1,:)~=0;
ch1 = ch1(2,thisTimebinHasPhotons);

% At this point we have kept only the timestamps with photons in them - and
% the number of photons is now no longer important.
timeBS = timeBS + toc;

% Apply detection: Call Transmitted detections channel 0 and Reflected
% detections channel 1
tic
ch0 = deadtime(ch0,deadTime); %take into account DT
ch0(2,:) = 0; %channel marker (as ch0-ch1 is positive delta tau)

ch1 = deadtime(ch1,deadTime); %take into account DT
ch1(2,:) = 1; %channel marker (as ch1-ch0 is negative delta tau)

timeDT = timeDT + toc;
% -------------------------------------------------------------------------
tic
% Calculate coincidences:
% If ith ch0 was detected before ch1 (i.e. ch0_timestamp < ch1_timestamp) the
% channel marker will be 0.5-(-0.5) = 1, thus interphotonDelays(2,i) =1.
% Likewise if ch1 was detected before ch0, will have -0.5-0.5 = -1.

interphotonDelays = diff(sortrows([ch0 ch1]',1),1,1)';
% Disregard interphoton delays between consecutive same channel detections:
sameChannelDelays = interphotonDelays(2,:)==0;
interphotonDelays(:,sameChannelDelays) = []; %this removes same channel coincidences

coincidencesFromThisLoop = interphotonDelays(1,:).*interphotonDelays(2,:);


coincidences = [coincidences, coincidencesFromThisLoop];
coincidencesOverTwiceRange = abs(coincidences) >(2*plottingRange);
coincidences(coincidencesOverTwiceRange) = [];

clear pulsetrain ch0 ch1 interphotonDelays coincidencesFromThisLoop
timeSortingCoinc = timeSortingCoinc + toc;
end

% -------------------------------------------------------------------------
tic
% Plot histogram of coincidences:
thisHistericalHistogram = figure;
lbx = 'Relative Delay [ns]';
lby = 'Coincidences';
lbtit = strcat('Two single photon emitters: \bf{\tau_1}=', num2str(tauDecay1),...
    'ns, \bf{\tau_2}=', num2str(tauDecay2));
h1 = histogram(coincidences, -plottingRange: bin :plottingRange);

h1.FaceColor = [0.5, 0.5, 0.5];
fontsize = 20;
xlabel(lbx);
ylabel(lby);
title(lbtit);
set(gca,'FontSize',12);
% set(gcf, 'Position', [100, 300, 900, 700]);
pbaspect([9 7 1]);
% display('plotting:')
toc
display(timeInit)
display(timeLoss)
display(timeBS)
display(timeDT)
display(timeSortingCoinc)

% addpath('./singlePhotonEmitter')
% filename = strcat('singlePhoton_tauDecay1_', num2str(tauDecay1), 'ns.mat');
filename = strcat('twoSPE_tauDecay1_', num2str(tauDecay1),...
    'ns_tauDecay2_', num2str(tauDecay2), 'ns.mat');
save(filename, 'coincidences')
% filenameHist = strcat('histogram_tauDecay1_', num2str(tauDecay1), 'ns.png');
filenameHist = strcat('histogramTwoSPE_tauDecay1_', num2str(tauDecay1),...
    'ns_tauDecay2_', num2str(tauDecay2), 'ns.png');
saveas(thisHistericalHistogram,filenameHist);

    end
end