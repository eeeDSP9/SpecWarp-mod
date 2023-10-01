myFolder = 'C:\Users\user\MATLAB\Spectral Warping\SWcode'
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end
filePattern = fullfile(myFolder, '**/*.wav'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
   file = theFiles(k).name
    [filepath,name,ext] = fileparts(file);
    [d,sr] = audioread(file);
    [a,g,e] = lpcfit(d,20);
    alpha = 0.05;
    [bhat, ahat]  = warppoles(a, alpha);

    dw = filter(bhat(1,:), 1, lpcsynth(ahat, g, e));

    soundsc(d,sr);
    soundsc(dw, sr);
    write_filename = name+"sw.wav";
    write_dir = 'C:\Users\user\MATLAB\Spectral Warping\sw(+0.05)\';
    audiowrite(write_dir+write_filename,dw,sr)
    clear sound;
end

    