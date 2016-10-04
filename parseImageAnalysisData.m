function out = parseImageAnalysisData(data, gene)
%% The Overview

% By Derek Fulton
% Last Modified: July 30 2016

% This program takes a data cell of data from Josh Lawrimore's
% image_analysis.m GUI and parses it into another data cell that can be
% easily used by MATLAB's classification learner for both easy data
% visualization and Machine Learning Techniques. 

% It is important to label the gene/mutation for the clustering category.

% For each 30s time point, this code calculates the following:

% Centroid x
% Centroid y
% Change in Centroid Position from Previous Time Point


% SPB1 x
% SPB1 y
% SPB1 z
% SPB2 x
% SPB2 y
% SPB2 z
% Spindle Length
% Change in Spindle Length from Previous Time Point's Spindle Length


% Major Axis Length (nm)
% Minor Axis Length (nm)
% Aspect Ratio
% Change in Major Axis Length (nm)
% Change in Minor Axis Length (nm)


% kMT length 1 
% kMT lenght 2
% Change in kMT 1 
% Change in kMT 2
% Difference length between kMT 1 and kMT 2

%% Example

% parsed_data = parseImageAnalysisData(data, 'sir2delete')
% where data is a data cell and <strain/mutation> is the unique identifier
% which classificationLearner uses to classify data into a mutation.

%% The Setup
centroid_column = 5; 
spb1_column = 7;
spb2_column = 8;
major_axis_column = 2;
minor_axis_column = 3;
kmt1_column = 9;
kmt2_column = 10;

rows = length(data);

goods = cell(rows, 11);
goods{1,1} = 'centroid_x';
goods{1,2} = 'centroid_y';
goods{1,3} = 'centroid_position_delta';
goods{1,4} = 'spb1_x';
goods{1,5} = 'spb1_y';
goods{1,6} = 'spb1_z';
goods{1,7} = 'spb2_x';
goods{1,8} = 'spb2_y';
goods{1,9} = 'spb2_z';
goods{1,10} = 'spindle_length';
goods{1,11} = 'spindle_length_delta';
goods{1,12} = 'major_axis';
goods{1,13} = 'minor_axis';
goods{1,14} = 'aspect_ratio';
goods{1,15} = 'major_axis_delta';
goods{1,16} = 'minor_axis_delta';
goods{1,17} = 'aspect_ratio_delta';
goods{1,18} = 'kmt1';
goods{1,19} = 'kmt2';
goods{1,20} = 'kmt_length_difference';
goods{1,21} = 'kmt1_stretch';
goods{1,22} = 'kmt2_stretch';
goods{1,23} = 'kmt_stretch_difference'; % Did one stretch more than the other?
goods{1,24} = 'Mutation';

%% The Loop
for i=2:rows
    %% The Centroids and Change in Centroid Position
    
    % Current Centroids
    if isempty(data{i, centroid_column})
        centroid_x = NaN;
        centroid_y = NaN;
    else
        centroid_x = data{i, centroid_column}(1); 
        centroid_y = data{i, centroid_column}(2); 
    end
    
    % Past Centroids 
    if isempty(data{i-1, centroid_column})
        centroid_x_past = NaN;
        centroid_y_past = NaN;
    else
        centroid_x_past = data{i-1, centroid_column}(1); 
        centroid_y_past = data{i-1, centroid_column}(2); 
    end
    
    % If either was empty or we're on the first row of the data, do nothing
    if isempty(centroid_x) || isempty(centroid_x_past) || i < 3
        centroid_position_delta = NaN;
    else % Calculate distance
        centroid_position_delta = sqrt( (centroid_x_past - centroid_x)^2 + (centroid_y_past - centroid_y)^2 );
    end
    
    
    %% Spindles, Spindle Length and Change in Spindle Length   
    
    % SPB 1 Present
    if isempty(data{i, spb1_column}) || (numel(data{i, spb1_column}) ~= 4)
        spb1_x = NaN; 
        spb1_y = NaN; 
        spb1_z = NaN; 
        spb1_intens = NaN;
    else
        spb1_x = data{i, spb1_column}(1); 
        spb1_y = data{i, spb1_column}(2); 
        spb1_z = data{i, spb1_column}(3); 
        spb1_intens = data{i, spb1_column}(4); 
    end
    
    % SPB 2 Present
    if (numel(data{i, spb2_column}) ~= 4)
        spb2_x = NaN; 
        spb2_y = NaN; 
        spb2_z = NaN; 
    else
        spb2_x = data{i, spb2_column}(1); 
        spb2_y = data{i, spb2_column}(2); 
        spb2_z = data{i, spb2_column}(3); 
    end
    
    % SPB 1 Past
    if isempty(data{i-1, spb1_column}) || (numel(data{i-1, spb1_column}) ~= 4)
        spb1_x_past = NaN; 
        spb1_y_past = NaN; 
        spb1_z_past = NaN; 
        spb1_intens_past = NaN;
    else
        spb1_x_past = data{i-1, spb1_column}(1); 
        spb1_y_past = data{i-1, spb1_column}(2); 
        spb1_z_past = data{i-1, spb1_column}(3); 
        spb1_intens_past = data{i-1, spb1_column}(4); 
    end
    
    % SPB 2 Past
    if isempty(data{i-1, spb2_column}) || (numel(data{i-1, spb2_column}) ~= 4)
        spb2_x_past = NaN; 
        spb2_y_past = NaN; 
        spb2_z_past = NaN; 
    else
        spb2_x_past = data{i-1, spb2_column}(1); 
        spb2_y_past = data{i-1, spb2_column}(2); 
        spb2_z_past = data{i-1, spb2_column}(3); 
    end
    
    
    % Spindle Length Calculations
    spindle_length = sqrt( (spb1_x - spb2_x)^2 + (spb1_y - spb2_y)^2 + (spb1_z - spb2_z)^2);
    spindle_length_past = sqrt( (spb1_x_past - spb2_x_past)^2 + (spb1_y_past - spb2_y_past)^2 + (spb1_z_past - spb2_z_past)^2);
    
    spindle_length_delta = spindle_length - spindle_length_past;
    
    %% Major Axis Length (Signal Length), Minor Axis Length, Aspect Ratio and Deltas
    
    % Current Major Axis
    if isempty(data{i, major_axis_column})
        major_axis = NaN;
    else
        major_axis = data{i, major_axis_column}(1); 
    end
    
    % Past Major Axis 
    if isempty(data{i-1, major_axis_column}) || (i < 3)
        major_axis_past = NaN;
    else
        major_axis_past = data{i-1, major_axis_column}(1); 
    end
    
    % Current Minor Axis
    if isempty(data{i, minor_axis_column})
        minor_axis = NaN;
    else
        minor_axis = data{i, minor_axis_column}(1); 
    end
    
    % Past Major Axis 
    if isempty(data{i-1, minor_axis_column}) || (i < 3)
        minor_axis_past = NaN;
    else
        minor_axis_past = data{i-1, minor_axis_column}(1); 
    end
    
    % Aspect Ratio
    if ~isempty(major_axis) && ~isempty(minor_axis)
        aspect_ratio = major_axis / minor_axis;
    else
        aspect_ratio = NaN;
    end
    
    % Past Aspect Ratio
    if isempty(major_axis_past) || isempty(minor_axis_past) || (i < 3)
        aspect_ratio_past = NaN;
    else
        aspect_ratio_past = major_axis_past / minor_axis_past;
    end
    
    % Major Axis (Signal Length) Delta
    if isempty(major_axis) || isempty(major_axis_past) || (i < 3)
        major_axis_delta = NaN;
    else
        major_axis_delta = major_axis - major_axis_past;
    end
    
    % Minor Axis Delta
    if isempty(minor_axis) || isempty(minor_axis_past) || (i < 3)
        minor_axis_delta = NaN;
    else
        minor_axis_delta = minor_axis - minor_axis_past;
    end
    
    % Aspect Ratio Delta
    if isempty(aspect_ratio) || isempty(aspect_ratio_past) || (i < 3)
        aspect_ratio_delta = NaN;
    else
        aspect_ratio_delta = aspect_ratio - aspect_ratio_past;
    end 
    
    %% The Kinetochore
    
    % 'Do they stretch evenly?'
    % kMT 1
    % past kMT 1
    % kMT 2
    % past kMT 2
    % Difference length between kMT 1 and kMT 2
    % Change in kMT 1 length (kmt1 stretch)
    % Change in kMT 2 length (kmt2 stretch)
    % Difference in stretch (are they stretching evenly?)
    
    % Current kMT 1 Length
    if isempty(data{i, kmt1_column})
        kmt1 = NaN;
    else
        kmt1 = data{i, kmt1_column}(1); 
    end
    
    % Past kMT 1 Length
    if isempty(data{i-1, kmt1_column}) || (i < 3)
        kmt1_past = NaN;
    else
        kmt1_past = data{i-1, kmt1_column}(1); 
    end
    
    % Current kMT 2 Length
    if isempty(data{i, kmt2_column})
        kmt2 = NaN;
    else
        kmt2 = data{i, kmt2_column}(1); 
    end
    
    % Past kMT 2 Length
    if isempty(data{i-1, kmt2_column}) || (i < 3)
        kmt2_past = NaN;
    else
        kmt2_past = data{i-1, kmt2_column}(1); 
    end
    
    % Difference Between Kinetochore Microtubule Lengths (Is one longer
    % than the other?)
    if isempty(kmt1) || isempty(kmt2)
        kmt_length_difference = NaN;
    else
        kmt_length_difference = abs( kmt1 - kmt2);
    end
    
    % Change in kMT 1 length (Stretch)
    if isempty(kmt1) || isempty(kmt1_past) || (i < 3)
        kmt1_stretch = NaN;
    else
        kmt1_stretch = kmt1 - kmt1_past;
    end
    
    % Change in kMT 2 length (Stretch)
    if isempty(kmt1) || isempty(kmt2_past) || (i < 3)
        kmt2_stretch = NaN;
    else
        kmt2_stretch = kmt2 - kmt2_past;
    end
    
    % Change in Difference Between Kinetochore Microtubule Length (Are they
    % stretching evenly?)
    if isempty(kmt1_stretch) || isempty(kmt2_stretch) || (i < 3)
       kmt_stretch_difference = NaN;
    else
       kmt_stretch_difference = kmt1_stretch - kmt2_stretch;
    end
    %% The Assignments
    goods{i,1} = {centroid_x};
    goods{i,2} = {centroid_y};
    goods{i,3} = {centroid_position_delta};
    goods{i,4} = {spb1_x};
    goods{i,5} = {spb1_y};
    goods{i,6} = {spb1_z};
    goods{i,7} = {spb2_x};
    goods{i,8} = {spb2_y};
    goods{i,9} = {spb2_z};
    goods{i,10} = {spindle_length};
    goods{i,11} = {spindle_length_delta};
    goods{i,12} = {major_axis};
    goods{i,13} = {minor_axis};
    goods{i,14} = {aspect_ratio};
    goods{i,15} = {major_axis_delta};
    goods{i,16} = {minor_axis_delta};
    goods{i,17} = {aspect_ratio_delta};
    goods{i,18} = {kmt1};
    goods{i,19} = {kmt2};
    goods{i,20} = {kmt_length_difference};
    goods{i,21} = {kmt1_stretch};
    goods{i,22} = {kmt2_stretch};
    goods{i,23} = {kmt_stretch_difference};
    goods{i,24} = {gene};

end
    
%% Parse the Table
vars = {'centroid_x', 'centroid_y', 'centroid_position_delta', 'spb1_x', 'spb1_y', 'spb1_z', 'spb2_x', 'spb2_y', 'spb2_z', 'spindle_length', 'spindle_length_delta', 'major_axis', 'minor_axis', 'aspect_ratio', 'major_axis_delta', 'minor_axis_delta', 'aspect_ratio_delta', 'kmt1', 'kmt2', 'kmt_length_difference', 'kmt1_stretch', 'kmt2_stretch', 'kmt_stretch_difference', 'Mutation'}

newcell = goods(2:end, :);

thetable = cell2table(newcell, 'VariableNames', vars);
writetable(thetable, 'temp.csv');
thetable = readtable('temp.csv'); % So we can concatenate them later
out = thetable;

end