function [] = LK_report_ft_timelockstatistics(outFT)
%
% LK_report_ft_timelockstatistics reports significance information
% contained in the output from ft_timelockstatistics.
%
% Use as: [] = LK_report_ft_timelockstatistics(outFT);
%
% Input: outFT, which is the output structure from ft_timelockstatistics.
%
% Lukas Kunz, 2022

% report positive cluster(s)
if isfield(outFT, 'posclusters') && ~isempty(outFT.posclusters)
    for iCluster = 1:max(outFT.posclusterslabelmat)
        fprintf('\nReporting fieldtrip output. Positive cluster #%d from %.3f to %.3f: t = %.3f, p = %.3f.\n', iCluster, ...
            min(outFT.time(outFT.posclusterslabelmat == iCluster)), max(outFT.time(outFT.posclusterslabelmat == iCluster)), ...
            outFT.posclusters(iCluster).clusterstat, outFT.posclusters(iCluster).prob);
    end
end

% report negative cluster(s)
if isfield(outFT, 'negclusters') && ~isempty(outFT.negclusters)
    for iCluster = 1:max(outFT.negclusterslabelmat)
        fprintf('\nReporting fieldtrip output. Negative cluster #%d from %.3f to %.3f: t = %.3f, p = %.3f.\n', iCluster, ...
            min(outFT.time(outFT.negclusterslabelmat == iCluster)), max(outFT.time(outFT.negclusterslabelmat == iCluster)), ...
            outFT.negclusters(iCluster).clusterstat, outFT.negclusters(iCluster).prob);
    end
end