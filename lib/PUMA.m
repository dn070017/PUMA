% 2017-06-16
% marieke kuijjer
% script by cho-yi
% OBS! i modified this to keep TFCoop fixed
% this will not run mirs and tfs together, but should be fine for running mirs alone
%
function RegNet = PANDA(RegNet, GeneCoReg, TFCoop, alpha)
    [NumTFs, NumGenes] = size(RegNet);
    disp('Learning Network!');
    tic;
    step = 0;
    hamming = 1;
    while hamming > 0.001
        R = Tfunction(TFCoop, RegNet);
        A = Tfunction(RegNet, GeneCoReg);
        W = (R + A) * 0.5; % was W = (R + A) * 0.5;
        hamming = mean(abs(RegNet(:) - W(:)));
        RegNet = (1 - alpha) * RegNet + alpha * W;

        if hamming > 0.001
            PPI = Tfunction(RegNet);
            PPI = UpdateDiagonal(PPI, NumTFs, alpha, step);            
            %%% OBS! below changed: only update diagonal
            TFCoop(eye(size(TFCoop))~=0)=(1 - alpha) * diag(TFCoop) + alpha * diag(PPI); % only update the diagonal

            CoReg2 = Tfunction(RegNet');
            CoReg2 = UpdateDiagonal(CoReg2, NumGenes, alpha, step);
            GeneCoReg = (1 - alpha) * GeneCoReg + alpha * CoReg2;
        end

        disp(['Step#', num2str(step), ', hamming=', num2str(hamming)]);
        step = step + 1;
        clear R A W CoReg2;  % release memory for next step
    end
    runtime = toc;
    fprintf('Running PANDA on %d Genes and %d TFs took %f seconds!\n', NumGenes, NumTFs, runtime);
end






