function out = LK_IndicateThousands(in)
%
% LK_IndicateThousands adds thousands separators to a 1x1 array.
% 
% Example: out = LK_IndicateThousands(1234567)
%
% Lukas Kunz, 2023.

% add commas to indicate thousands
import java.text.*
v = DecimalFormat;
out = char(v.format(in));
out = strrep(out, '.', ','); % convert dots into commas