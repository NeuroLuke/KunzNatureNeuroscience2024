function chans = LK_GetChans(chanDir)
%
% LK_GetChans gets all channels within a specific directory, sorted in
% ascending order.
%
% Use as: chans = LK_GetChans(chanDir);
%
% See also dir, replace, lettersPattern
%
% Lukas Kunz, 2021

% all channels in the directory
chans   = dir(chanDir);

% sort all channels in ascending order
[~, I]  = sort(cellfun(@str2double, replace({chans.name}, lettersPattern, '')));
chans   = chans(I);