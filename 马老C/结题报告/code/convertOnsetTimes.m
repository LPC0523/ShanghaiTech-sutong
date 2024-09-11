% %---------------------%
% % Convert Onset Times %
% %---------------------%
% 
% % Converts timing files from BIDS format into a two-column format that can
% % be read by SPM
% 
% % The columns are:
% % 1. Onset (in seconds); and
% % 2. Duration (in seconds
% 
% 
% % Run this script from the directory that contains all of your subjects
% % (i.e., the Flanker directory)

subjects = [01 03 04 08 09 11 12 13 15 17 18 19 20 21 22 23 26 28 30]; % Replace with a list of all of the subjects you wish to analyze


for subject=subjects
    
    subject = num2str(subject, '%02d'); % Zero-pads each number so that the subject ID is 2 characters long

    cd(['sub-' subject '\func']) % Navigate to the subject's directory

    Run1_onsetTimes = tdfread(['sub-' subject '_task-cpshop_run-01_events.tsv'], '\t'); % Read onset times file
    Run1_onsetTimes.trial_type = string(Run1_onsetTimes.trial_type); % Convert char array to string array, to make logical comparisons easier

    Run1_together = [];
    Run1_alone = [];

    for i = 1:length(Run1_onsetTimes.onset)
        if strtrim(Run1_onsetTimes.trial_type(i,:)) == 'together'
            Run1_together = [Run1_together; Run1_onsetTimes.onset(i,:) Run1_onsetTimes.duration(i,:)];
        elseif strtrim(Run1_onsetTimes.trial_type(i,:)) == 'alone'
            Run1_alone = [Run1_alone; Run1_onsetTimes.onset(i,:) Run1_onsetTimes.duration(i,:)];
        end
    end

    Run2_onsetTimes = tdfread(['sub-' subject '_task-cpshop_run-02_events.tsv'], '\t'); % Read onset times file
    Run2_onsetTimes.trial_type = string(Run2_onsetTimes.trial_type); % Convert char array to string array, to make logical comparisons easier

    Run2_together = [];
    Run2_alone = [];

    for i = 1:length(Run2_onsetTimes.onset)
        if strtrim(Run2_onsetTimes.trial_type(i,:)) == 'together'
            Run2_together = [Run2_together; Run2_onsetTimes.onset(i,:) Run2_onsetTimes.duration(i,:)];
        elseif strtrim(Run2_onsetTimes.trial_type(i,:)) == 'alone'
            Run2_alone = [Run2_alone; Run2_onsetTimes.onset(i,:) Run2_onsetTimes.duration(i,:)];
        end
    end
    
    Run3_onsetTimes = tdfread(['sub-' subject '_task-cpshop_run-03_events.tsv'], '\t'); % Read onset times file
    Run3_onsetTimes.trial_type = string(Run3_onsetTimes.trial_type); % Convert char array to string array, to make logical comparisons easier

    Run3_together = [];
    Run3_alone = [];

    for i = 1:length(Run3_onsetTimes.onset)
        if strtrim(Run3_onsetTimes.trial_type(i,:)) == 'together'
            Run3_together = [Run3_together; Run3_onsetTimes.onset(i,:) Run3_onsetTimes.duration(i,:)];
        elseif strtrim(Run3_onsetTimes.trial_type(i,:)) == 'alone'
            Run3_alone = [Run3_alone; Run3_onsetTimes.onset(i,:) Run3_onsetTimes.duration(i,:)];
        end
    end

    Run4_onsetTimes = tdfread(['sub-' subject '_task-cpshop_run-04_events.tsv'], '\t'); % Read onset times file
    Run4_onsetTimes.trial_type = string(Run4_onsetTimes.trial_type); % Convert char array to string array, to make logical comparisons easier

    Run4_together = [];
    Run4_alone = [];

    for i = 1:length(Run4_onsetTimes.onset)
        if strtrim(Run4_onsetTimes.trial_type(i,:)) == 'together'
            Run4_together = [Run4_together; Run4_onsetTimes.onset(i,:) Run4_onsetTimes.duration(i,:)];
        elseif strtrim(Run4_onsetTimes.trial_type(i,:)) == 'alone'
            Run4_alone = [Run4_alone; Run4_onsetTimes.onset(i,:) Run4_onsetTimes.duration(i,:)];
        end
    end

    % Save timing files into text files

    save('together_run1.txt', 'Run1_together', '-ASCII');
    save('together_run2.txt', 'Run2_together', '-ASCII');
    save('together_run3.txt', 'Run3_together', '-ASCII');
    save('together_run4.txt', 'Run4_together', '-ASCII');


    save('alone_run1.txt', 'Run1_alone', '-ASCII');
    save('alone_run2.txt', 'Run2_alone', '-ASCII');
    save('alone_run3.txt', 'Run3_alone', '-ASCII');
    save('alone_run4.txt', 'Run4_alone', '-ASCII');

    % Go back to Flanker directory

    cd ../..

end
