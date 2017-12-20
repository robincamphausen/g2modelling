function [ output ] = deadtime( input, deadTime )
%deadtime function: Removes timestamps that would not be detected due to
%detector lifetime
% Author: Robin Camphausen
% Last modified: 12/12/2017

% -Very first timestamp is kept.
% -Calculates the interphoton delays, i.e. the diff between consecutive
% detections.
% -Calculates the time since the last detection: this initialises at zero
% and then while looping through the diff array always adds on the latest
% diff
% -If the time since last detection is smaller than deadtime, discard that
% timestamp, if larger keep the timestamp and reset time
% -------------------------------------------------------------------------

timeSinceLastDetection = 0;
interphotonDelays = abs(diff(input));

for i = 1:length(interphotonDelays)
    timeSinceLastDetection = timeSinceLastDetection + interphotonDelays(i);
    if timeSinceLastDetection < deadTime
        input(i+1) = -1; %mark this timestamp for deletion
    elseif timeSinceLastDetection >= deadTime
        timeSinceLastDetection = 0; %photon detected and deadtime reset
    end
end

keepThese = input>=0;
output = input(keepThese);
end