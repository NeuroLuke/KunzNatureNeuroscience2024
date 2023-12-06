function [] = LK_save(saveFile, saveVariable, saveFlag)
%
% LK_save replaces the MATLAB function save to be able to save variables
% within a parfor loop.
%
% Input:
%   saveFile        --> path and name under which the variable shall be
%                       saved
%   saveVariable    --> variable to save (without quotation marks/ticks)
%   saveFlag        --> whether to use the '-v7.3' option when saving large
%                       variables.
%
% Lukas Kunz, 2021

% original variable name of the second input
varName             = inputname(2);

% use a structure to save the desired variable, which you add as a field to
% the structure
saveVar             = [];
saveVar.(varName)   = saveVariable; % assign the actual data to the field with the desired name of the variable to save

% save the fields of the structure "saveVar" as individual variable in the
% file "saveFile"
if nargin > 2
    save(saveFile, '-struct', 'saveVar', saveFlag);
else
    save(saveFile, '-struct', 'saveVar');
end
