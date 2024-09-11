% List of open inputs
nrun = 19; % enter the number of runs here
jobfile = {'D:\shopping interaction\estimate_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
