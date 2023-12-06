function number = LK_letter2num(letter)
%
% LK_num2letter converts one or more alphabetical letters (e.g., A) into
% their corresponding number (e.g., 1). The output number is given as
% int16. It is not possible to have both uppercase and lowercase letters in
% the input letter at the same time.
%
% Examples: a --> 1; A --> 1; bb --> [2, 2]; ZZZ --> [26, 26, 26].
%
% Use as: number = LK_letter2num(letter)
%
% Lukas Kunz, 2022

% check whether all letters are either uppercase or lowercase
if all(int16(letter) > 90)
    bUppercase = false;
elseif all(int16(letter) <= 90)
    bUppercase = true;
else
    bUppercase = nan;
end

% convert letters into numbers
if isnan(bUppercase)
    error('Letter contains both uppercase and lowercase letters.');
elseif bUppercase == false
    number = int16(letter) - int16('a') + 1;
elseif bUppercase == true
    number = int16(letter) - int16('A') + 1;
end