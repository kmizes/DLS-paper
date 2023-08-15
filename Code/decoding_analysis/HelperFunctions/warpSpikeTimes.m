function [plotspikes, wpsth] = warpSpikeTimes(sptimes, events, templ, histedges, samprate)

% Loop across trials in this session and plot spike times warped to template
ntrials = size(events,1);
nevents = size(events,2);
devents = diff(events,1,2);
if nargin < 5
    samprate = 3e4;
end
sptimes = double(sptimes);

wpsth = zeros(ntrials, length(histedges));
plotspikes = [];

for t = 1 : ntrials
    spikes = sptimes((sptimes >= (events(t,1)+histedges(1)*samprate)) & (sptimes <= (events(t,1)+histedges(end)*samprate)));
    
    wspikes = zeros(size(spikes));
    for sp = 1 : length(spikes)
        pevent = find(events(t,:) <= spikes(sp), 1, 'last');
        if ~isempty(pevent) && pevent ~= nevents
            if pevent ~= 1
                toffset = sum(templ(1:pevent-1));
            else
                toffset = 0;
            end
            wspikes(sp) = toffset + ((spikes(sp)-events(t,pevent))*templ(pevent)/devents(t,pevent));
        elseif isempty(pevent)
            wspikes(sp) = spikes(sp)-events(t,1);
        else
            wspikes(sp) = spikes(sp)-(events(t,end)-sum(templ(1:end)));
        end
    end

    % Convert to seconds timescale
    wspikes = wspikes/samprate;
    wpsth(t, :) = histc(wspikes, histedges);

    plotspikes = [plotspikes [wspikes'; wspikes'*0+t]];
end

wpsth = wpsth(:,1:end-1);