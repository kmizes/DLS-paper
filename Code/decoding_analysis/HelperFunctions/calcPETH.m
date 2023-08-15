function upeth = calcPETH(sptimes, events, histedges, samprate, mem_lim)
% Arrange spikes into PSTH based on eventTimes and histedges. Samprate
% input optional (3e4 Hz by default).

if nargin <4
    samprate = 3e4;
end

if nargin < 5
    mem_lim = 15e9;
end

ntrials = size(events,1);
sptimes = double(sptimes);

upeth = zeros(ntrials, length(histedges));

if ntrials*length(sptimes)*2*8 > mem_lim
    for t = 1 : ntrials
        spikes = sptimes((sptimes >= (events(t,1)+histedges(1)*samprate)) & (sptimes <= (events(t,1)+histedges(end)*samprate)));
        spikes = (spikes - events(t,1))/samprate;
        upeth(t,:) = histc(spikes, histedges);
    end
    upeth = upeth(:,1:end-1);
else
    spikes = sptimes(sptimes >= min(events(:))+histedges(1)*samprate & sptimes <= (max(events(:))+histedges(end)*samprate));
    d_spev = (repmat(spikes, 1, ntrials) - repmat(events(:,1)', length(spikes), 1)) / samprate;
    [upeth, ~] = histc(d_spev, histedges);
    upeth = upeth(1:end-1,:)';
    
    if isempty(upeth)
        upeth = zeros(ntrials, length(histedges)-1);
    end    
end
    
