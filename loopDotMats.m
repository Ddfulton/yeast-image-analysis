function out = loopDotMats(strain)
%% Description
%   Loops through all *.mat files in the current directory and concatenates 
%   them into one massive table suitable for data mining

files = dir('*.mat');
files = {files.name};


for i=1:length(files)
    z = cellstr(files(i));
    z = z{1};
    load(z);
    title = strcat(strain, int2str(i));
    % title = parseImageAnalysisData(data_cell, 'sir21');
    expression = strcat(title, '=parseImageAnalysisData(data_cell,''', strain, ''')');
    eval(expression);
    
end

out = title;

clear z;
end