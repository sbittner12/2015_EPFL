signal_id = 'pair2';
pmin = 40;
pmax = pmin;
windowsize = 512;

% Load pascal result
file_path = 'E:\Sean\Code\iCohForFree\EEG_EMG_pair\';
load([file_path, signal_id, '.mat']);
load([file_path, signal_id, '_phase_randomized.mat']);

fname_original = [file_path, sprintf('%s-iCOH.txt', signal_id)];
fname_pr =  [file_path, sprintf('%s_phase_randomized-iCOH.txt', signal_id)];
icoh_normal_pascal = dlmread(fname_original, ' ');
icoh1_normal_pascal = squeeze(icoh_normal_pascal(1:2:512,5));
icoh2_normal_pascal = squeeze(icoh_normal_pascal(1:2:512,7));

icoh_pr_pascal = dlmread(fname_pr, ' ');
icoh1_pr_pascal = squeeze(icoh_pr_pascal(1:2:512,5));
icoh2_pr_pascal = squeeze(icoh_pr_pascal(1:2:512,7));

% Load matlab result
[bb_normal, dim, acmat] = iCOH(x_normal, pmin, windowsize);
[bb_pr, dim, acmat] = iCOH(x_pr, pmin, windowsize);
icoh1_normal_matlab = bb_normal(:,2);
icoh2_normal_matlab = bb_normal(:,3);
icoh1_pr_matlab = bb_pr(:,2);
icoh2_pr_matlab = bb_pr(:,3);



% Load ARfit incorporated result

[w, A_normal, C_normal, sbc, fpe, th] = arfit(x_normal, pmin, pmax);
[w, A_pr, C_pr, sbc, fpe, th] = arfit(x_pr, pmin, pmax);

acfunc_normal = getacfunc(C_normal, A_normal);
acfunc_pr = getacfunc(C_pr, A_pr);

A_w_normal = Autoco2FFT(acfunc_normal, windowsize);
A_w_pr = Autoco2FFT(acfunc_pr, windowsize);

icoh_normal_arfit = computeiCOH(C_normal, A_w_normal, windowsize);
icoh1_normal_arfit = squeeze(icoh_normal_arfit(2,1,:));
icoh2_normal_arfit = squeeze(icoh_normal_arfit(1,2,:));

icoh_pr_arfit = computeiCOH(C_pr, A_w_pr, windowsize);
icoh1_pr_arfit = squeeze(icoh_pr_arfit(2,1,:));
icoh2_pr_arfit = squeeze(icoh_pr_arfit(1,2,:));

figure;
subplot(2,3,1), plot(icoh1_normal_pascal);
title('Pascal EEG -> EMG');
subplot(2,3,4), plot(icoh2_normal_pascal);
title('Pascal EMG -> EEG');
subplot(2,3,2), plot(icoh1_normal_matlab);
title('Matlab EMG -> EEG');
subplot(2,3,5), plot(icoh2_normal_matlab);
title('Matlab EMG -> EEG');
subplot(2,3,3), plot(icoh1_normal_arfit);
title('Arfit EEG -> EMG');
subplot(2,3,6), plot(icoh2_normal_arfit);
title('Arfit EMG -> EEG');

figure;
subplot(2,3,1), plot(icoh1_pr_pascal);
title('Pascal EEG -> EMG');
subplot(2,3,4), plot(icoh2_pr_pascal);
title('Pascal EMG -> EEG');
subplot(2,3,2), plot(icoh1_pr_matlab);
title('Matlab EEG -> EMG');
subplot(2,3,5), plot(icoh2_pr_matlab);
title('Matlab EMG -> EEG');
subplot(2,3,3), plot(icoh1_pr_arfit);
title('Arfit EEG -> EMG');
subplot(2,3,6), plot(icoh2_pr_arfit);
title('Arfit EMG -> EEG');

