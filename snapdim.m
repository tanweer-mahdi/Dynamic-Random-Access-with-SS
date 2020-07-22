function [snaps,N] = snapdim(Y,m)

%% This function computes snapshots from signal observation for a given length
% *** Inputs ***
% Y = Observation vector
% m = snapshot length
% *** Outputs ***
% snaps = Snapshot matrix
% N = Number of snapshots

init = (1:m)';
lel = init;
while(1)
    lel = lel+1;
    temp = lel(:,end);
    if temp(end)>length(Y)
        break;
    end
    init = [init lel];
end

snaps = [];
for i=1:size(init,2)
    snaps = [snaps Y(init(:,i)) conj(flipud(Y(init(:,i))))];
end
end

        