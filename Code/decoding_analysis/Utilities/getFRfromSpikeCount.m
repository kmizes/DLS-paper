function fr = getFRfromSpikeCount(spikes, binsize, window)

% Estimates smoothed firing rate by integrating spike counts over some
% binsize (in ms!)

if nargin < 3
    window = 100; % ms!
end
fr = [];

T = ceil(window / binsize);
% pad
spikes = [spikes(1)*ones(1,round(T/2)) spikes spikes(end)*ones(1,round(T/2))];
% sliding window
for t = (round(T/2)+1):(length(spikes)-T+round(T/2)) % used to be -1 here?
    fr(t-round(T/2)) = sum(spikes((t-round(T/2)):(t+round(T/2)-1))) ./ (binsize*T) * 1000;
end

% smooth by window size?
fr = smooth(fr, T);

% readjust to full size of spike input?

%*** optional, try convolve with gaussian to get FR
% build kernel

% sz = binsize * 5; % 
% sigma = binsize * .5;
% x = linspace(-sz / 2, sz / 2, sz);
% gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
% gaussFilter = gaussFilter / sum (gaussFilter); % normalize

%figure;plot(conv(sptemp,gaussFilter,'same'))
%hold on;plot(find(sptemp),0.025*ones(1,length(find(sptemp))),'.')

end