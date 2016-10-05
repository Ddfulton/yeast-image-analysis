function varargout = population_GUI_v1_2(varargin)
% POPULATION_GUI_V1_2 MATLAB code for population_GUI_v1_2.fig
%      POPULATION_GUI_V1_2, by itself, creates a new POPULATION_GUI_V1_2 or raises the existing
%      singleton*.
%
%      H = POPULATION_GUI_V1_2 returns the handle to a new POPULATION_GUI_V1_2 or the handle to
%      the existing singleton*.
%
%      POPULATION_GUI_V1_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POPULATION_GUI_V1_2.M with the given input arguments.
%
%      POPULATION_GUI_V1_2('Property','Value',...) creates a new POPULATION_GUI_V1_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before population_GUI_v1_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to population_GUI_v1_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help population_GUI_v1_2

% Last Modified by GUIDE v2.5 08-Jul-2016 14:10:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @population_GUI_v1_2_OpeningFcn, ...
    'gui_OutputFcn',  @population_GUI_v1_2_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before population_GUI_v1_2 is made visible.
function population_GUI_v1_2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to population_GUI_v1_2 (see VARARGIN)

% Choose default command line output for population_GUI_v1_2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes population_GUI_v1_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = population_GUI_v1_2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function plane_slider_Callback(hObject, eventdata, handles)
% hObject    handle to plane_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% set slider_value to the plane_slider 'Value' GFP channel

%load in the existing centroids
centroid_points = handles.centroid_points;
% update axes3 first so that GFP channel is "selected" for thresholding
slider_value = ceil(get(handles.plane_slider,'Value'));
if isempty(centroid_points) == 0
    subplot(handles.axes3);
    hold on;
    imshow(handles.im_mat_2(:,:,slider_value),[]);
    plot(centroid_points(:,1), centroid_points(:,2), 'go');
    hold off;
 else
    subplot(handles.axes3);
    imshow(handles.im_mat_2(:,:,slider_value),[]);
end

if isempty(centroid_points) == 0
    subplot(handles.axes1);
    hold on;
    imshow(handles.im_mat(:,:,slider_value),[]);
    plot(centroid_points(:,1), centroid_points(:,2), 'go');
    hold off;
else
    subplot(handles.axes1);
    imshow(handles.im_mat(:,:,slider_value),[]);
end
%set slider text
plane_num = handles.plane_num;
set(handles.slider_text, 'String', strcat(num2str(...
    ceil(get(handles.plane_slider,'Value'))),...
    '/', num2str(plane_num)));
%set Stack text
stack_total = plane_num/str2double(get(handles.step_number,'String'));
current_stack = ceil(slider_value/...
    str2double(get(handles.step_number,'String')));
set(handles.stack_text,'String',strcat(num2str(current_stack),...
    num2str('/'), num2str(stack_total)));





% --- Executes during object creation, after setting all properties.
function plane_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plane_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function open_image_Callback(hObject, eventdata, handles)
% hObject    handle to open_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set working directory to current directory
wdir = pwd;

%select the file
[filename pathname] = uigetfile('*.tif', 'Choose the GFP image stack');
%read in the file using tif3Dread.m
im_mat = tif3Dread(strcat(pathname,filename));
%change image to double
im_mat = double(im_mat);
%set up data to axes1 plot
subplot(handles.axes1);
%show first frame of image stack
imshow(im_mat(:,:,1),[]);
%read in the number of frames
[~,~,plane_num] = size(im_mat);
handles.plane_num = plane_num;
%set slider
slider_step = 1/plane_num;
set(handles.plane_slider, 'SliderStep', [slider_step slider_step],...
    'Min', 1, 'Max', plane_num, 'Value', 1);
handles.im_mat = im_mat;
%set slider text
set(handles.slider_text, 'String', strcat(num2str(...
    get(handles.plane_slider,'Value')),...
    '/', num2str(plane_num)));
%set Stack text
stack_total = plane_num/str2double(get(handles.step_number,'String'));
current_stack = 1;
set(handles.stack_text,'String',strcat(num2str(current_stack),...
    num2str('/'), num2str(stack_total)));

%open the RFP stack in axes3
[filename pathname] = uigetfile('*.tif', 'Choose the RFP image stack');
im_mat_2 = tif3Dread(strcat(pathname,filename));
im_mat_2 = double(im_mat_2);
handles.im_mat_2 = im_mat_2;
subplot(handles.axes3);
imshow(im_mat_2(:,:, ceil(get(handles.plane_slider,'Value'))),[]);
%set up the cell array that will hold the image stack data
data_cell = {'Plane', 'Major Axis (nm)','Minor Axis (nm)', 'Aspect Ratio',...
    'Centroid', '3D Spindle (nm)','SPB 1', 'SPB 2', 'kMT 1', 'kMT2',...
    'Plasmid Split?','Plasmid Image'};
handles.data_cell = data_cell;
record_num = 1;
handles.record_num = record_num;
%instantiate the handles.centroid_points
handles.centroid_points = [];

% update the handles data structure
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in FG_button.
function FG_button_Callback(hObject, eventdata, handles)
% hObject    handle to FG_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = imrect;
pos_fg = wait(h);
pos_fg = ceil(pos_fg);
handles.pos_fg = pos_fg;
im_mat_fg = handles.im_mat(pos_fg(2):pos_fg(2)+(pos_fg(4)),pos_fg(1):...
    pos_fg(1)+pos_fg(3), ceil(get(handles.plane_slider,'Value')));
subplot(handles.axes2);
imshow(im_mat_fg, []);
mean(im_mat_fg(:));
k = waitforbuttonpress;

%subtract out the average BG from FG
im_mat_fg = im_mat_fg - handles.avg_bg;

%convert the negative numbers to zeros
im_mat_fg(im_mat_fg<0) = 0;
subplot(handles.axes2);
imshow(im_mat_fg,[]);
k = waitforbuttonpress;

%create a binary based on multithresh
im_thresh = multithresh(im_mat_fg);
im_bin = im_mat_fg > im_thresh;
subplot(handles.axes2);
imshow(im_bin, []);
handles.im_bin = im_bin;


%pad the binary image back to the orginal size
[im_y, im_x, ~] = size(handles.im_mat);
x_pad_1 = zeros(pos_fg(4) + 1, pos_fg(1) - 1);
im_bin_resize = horzcat(x_pad_1,im_bin);
[~,im_bin_resize_width] = size(im_bin_resize);
x_pad_2 = zeros(pos_fg(4) +1, im_x - im_bin_resize_width);
im_bin_resize = horzcat(im_bin_resize, x_pad_2);

y_pad_1 = zeros(pos_fg(2) - 1, im_x);
im_bin_resize = vertcat(y_pad_1,im_bin_resize);
[im_bin_resize_height, ~] = size(im_bin_resize);
y_pad_2 = zeros(im_y - im_bin_resize_height, im_x);
im_bin_resize = vertcat(im_bin_resize, y_pad_2);
slider_value = ceil(get(handles.plane_slider,'Value'));
handles.plasmid_mask = im_bin_resize;
handles.plasmid_plane = slider_value;

% Save the plasmid image, Note: this does not have the bg sub
size(handles.im_mat(:,:,slider_value))
size(im_bin_resize)
im_check = handles.im_mat(:,:,slider_value) .* im_bin_resize;
handles.im_check = im_check;


%get the centroid of the resized binary image
if get(handles.no_plasmid_check,'Value') == 0
    centroid_struct = regionprops(im_bin_resize,'centroid');
    centroid = centroid_struct(1);
    handles.centroid = struct2array(centroid);
    %plot the centroid
    subplot(handles.axes1);
    hold on
    plot(handles.centroid(1,1), handles.centroid(1,2), 'go');
    hold off
    subplot(handles.axes3);
    hold on
    plot(handles.centroid(1,1), handles.centroid(1,2), 'go');
    hold off
else
    handles.centroid = [];
    handles.plasmid_plane = [];
end

handles.centroid_points = [handles.centroid_points; struct2array(centroid)];
guidata(hObject,handles);



% --- Executes on button press in sub_bg.
function sub_bg_Callback(hObject, eventdata, handles)
% hObject    handle to sub_bg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.im_mat(:,:,ceil(get(handles.plane_slider,'Value'))) = ...
    handles.im_mat(:,:,ceil(get(handles.plane_slider,'Value'))) - ...
    handles.avg_bg;
max(max(handles.im_mat(:,:,ceil(get(handles.plane_slider,'Value')))))
%set up data to axes1 plot
subplot(handles.axes1);
imshow(handles.im_mat(:,:,ceil(get(handles.plane_slider,'Value'))),[]);
set(handles.bg_text, 'String', 'Applied BG sub');


% --- Executes on button press in major_axis.
function major_axis_Callback(hObject, eventdata, handles)
% hObject    handle to major_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get major axis length
bw_props = regionprops(handles.im_bin,'MajorAxisLength');
pixel_size = str2double(get(handles.pixel_size, 'String'));
axis_length = bw_props.MajorAxisLength * pixel_size;
%get minor axis length
bw_props_2 = regionprops(handles.im_bin,'MinorAxisLength');
minor_axis = bw_props_2.MinorAxisLength * pixel_size;
%calculate aspect ratio
aspect_rat = axis_length/minor_axis;
%puth major and minor axis to handles structure
handles.axis_length = axis_length;
handles.minor_axis = minor_axis;
handles.aspect_rat = aspect_rat;
guidata(hObject,handles);
%set the text boxes
set(handles.axis_tag, 'String', num2str(round(axis_length)));
set(handles.minor_tag,'String', num2str(round(minor_axis)));
set(handles.aspect_rat_tag,'String', num2str(aspect_rat));



% --- Executes on button press in spb_1_button.
function spb_1_button_Callback(hObject, eventdata, handles)
% hObject    handle to spb_1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Switch to axes3
subplot(handles.axes3);
%wait for user button press on SPB
k = waitforbuttonpress;
click_point = round(get(gca,'CurrentPoint'));
hold on
plot(click_point(1,1),click_point(1,2), '-rx');
hold off
%get the current plane to find Z position
current_plane=ceil(get(handles.plane_slider, 'Value'));
%figure out where we are in the z-stack
%get the step number information from the edit text box
step_number = str2double(get(handles.step_number,'String'));
%modulo of current step by step number
stack_mod = mod(ceil(get(handles.plane_slider,'Value')), step_number);
if stack_mod == 0
    stack_min = (1-step_number) + current_plane;
    stack_max = current_plane;
else
    stack_min = (1-stack_mod) + current_plane;
    stack_max = abs(stack_mod - step_number) + current_plane;
end

% let me display the search area and the resulting pixel
% on the command line to check
% current_plane
% stack_min
% stack_max

% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
if click_point(1,1)<5
    click_point(1,1)=5;
end
if click_point(1,2)<5
    click_point(1,2)=5;
end

%correcting if click region is close to maximum
im_mat = handles.im_mat;
[y,x,~] = size(im_mat);
if click_point(1,1) > (x-4)
    click_point(1,1) = x-4;
end
if click_point(1,2) > (y-4)
    click_point(1,2) = y-4;
end
crop_im = handles.im_mat_2((click_point(1,2)-4):(click_point(1,2)+4),...
    (click_point(1,1)-4):(click_point(1,1)+4),stack_min:stack_max);
max_pix = max(crop_im(:));
[max_y, max_x, max_z] = ind2sub(size(crop_im),find(crop_im==max_pix));
if length(max_x)~= 1
    max_x=max_x(1);
    max_y=max_y(1);
    max_z=max_z(1);
end
%convert from crop_im coordinates back to the im_mat_2 coordinates
max_x_im = max_x + click_point(1,1)-4 - 1;
max_y_im = max_y + click_point(1,2)-4 - 1;
max_z_im = max_z + stack_min - 1;
test_max = handles.im_mat_2(max_y_im, max_x_im, max_z_im);
% %check that the intensity values match by displaying the crop and
% %the original
% % test_max
% % max_pix
% save the X Y and Plane data for SPB1 to the handles structure
spb1 = [max_x_im, max_y_im, max_z_im, max_pix];
handles.spb1 = spb1;
guidata(hObject,handles);

% --- Executes on button press in spb_2_button.
function spb_2_button_Callback(hObject, eventdata, handles)
% hObject    handle to spb_2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Switch to axes3
subplot(handles.axes3);
%wait for user button press on SPB
k = waitforbuttonpress;
click_point = round(get(gca,'CurrentPoint'));
hold on
plot(click_point(1,1),click_point(1,2), '-rx');
hold off
%get the current plane to find Z position
current_plane=ceil(get(handles.plane_slider, 'Value'));
%figure out where we are in the z-stack
%get the step number information from the edit text box
step_number = str2double(get(handles.step_number,'String'));
%modulo of current step by step number
stack_mod = mod(ceil(get(handles.plane_slider,'Value')), step_number);
if stack_mod == 0
    stack_min = (1-step_number) + current_plane;
    stack_max = current_plane;
else
    stack_min = (1-stack_mod) + current_plane;
    stack_max = abs(stack_mod - step_number) + current_plane;
end

% let me display the search area and the resulting pixel
% on the command line to check
% current_plane
% stack_min
% stack_max
% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
if click_point(1,1)<5
    click_point(1,1)=5;
end
if click_point(1,2)<5
    click_point(1,2)=5;
end

%correcting if click region is close to maximum
im_mat = handles.im_mat;
[y,x,~] = size(im_mat);
if click_point(1,1) > (x-4)
    click_point(1,1) = x-4;
end
if click_point(1,2) > (y-4)
    click_point(1,2) = y-4;
end

% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
crop_im = handles.im_mat_2((click_point(1,2)-4):(click_point(1,2)+4),...
    (click_point(1,1)-4):(click_point(1,1)+4),stack_min:stack_max);
max_pix = max(crop_im(:));
[max_y, max_x, max_z] = ind2sub(size(crop_im),find(crop_im==max_pix));
if length(max_x)~= 1
    max_x=max_x(1);
    max_y=max_y(1);
    max_z=max_z(1);
end
%convert from crop_im coordinates back to the im_mat_2 coordinates
max_x_im = max_x + click_point(1,1)-4 - 1;
max_y_im = max_y + click_point(1,2)-4 - 1;
max_z_im = max_z + stack_min - 1;
test_max = handles.im_mat_2(max_y_im, max_x_im, max_z_im);
% %check that the intensity values match by displaying the crop and
% %the original
% test_max
% max_pix
spb2 = [max_x_im, max_y_im, max_z_im, max_pix];
handles.spb2 = spb2;
guidata(hObject,handles);



function step_size_Callback(hObject, eventdata, handles)
% hObject    handle to step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_size as text
%        str2double(get(hObject,'String')) returns contents of step_size as a double


% --- Executes during object creation, after setting all properties.
function step_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_number_Callback(hObject, eventdata, handles)
% hObject    handle to step_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_number as text
%        str2double(get(hObject,'String')) returns contents of step_number as a double


% --- Executes during object creation, after setting all properties.
function step_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixel_size_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_size as text
%        str2double(get(hObject,'String')) returns contents of pixel_size as a double


% --- Executes during object creation, after setting all properties.
function pixel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calc_button.
function calc_button_Callback(hObject, eventdata, handles)
% hObject    handle to calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%make contigent on plasmid detected or not

%load in the SPB positions
if get(handles.no_spb1_check, 'Value') == 0
    spb1 = handles.spb1;
else
    spb1 = [];
    handles.spb1 = spb1;
end

if get(handles.no_spb2_check, 'Value') == 0
    spb2 = handles.spb2;
    handles.spb2 = spb2;
else
    spb2 = [];
    handles.spb2 = spb2;
end
%load in the step size and pixel size
pixel_size = str2double(get(handles.pixel_size, 'String'));
step_size = str2double(get(handles.step_size, 'String'));
if get(handles.no_spb1_check,'Value') == 1 || ...
        get(handles.no_spb2_check, 'Value') == 1 || ...
        get(handles.split_check, 'Value') == 1
    %set the spindle_3D value to blank
    spindle_3D = [];
    handles.spindle_3D = spindle_3D;
    % Display that the spindle length cannot be calculated
    set(handles.spindle_nm, 'String', 'NA');
else
    %calculate the 3D spindle length
    spindle_3D = sqrt(((spb1(1)-spb2(1))*pixel_size)^2 +...
        (((spb1(2) - spb2(2))*pixel_size)^2) +...
        (((spb1(3)-spb2(3))*step_size)^2));
    handles.spindle_3D = spindle_3D;
    % Display the spindle length in the GUI
    set(handles.spindle_nm, 'String', num2str(round(spindle_3D)));
end

%loop through all the 1 positions in the binary and calculate
%the distance from SPB1 to each pixel
%load in the resized plasmid binary image
im_bin_resize = handles.plasmid_mask;
%load in the plasmid plane
plasmid_plane = handles.plasmid_plane;

bin_ind = find(im_bin_resize);
[y, x] = ind2sub(size(im_bin_resize),bin_ind);
%merge x and y plasmid coordinates with plane information
plane_array = ones(length(bin_ind),1) * plasmid_plane;
plasmid_coords = horzcat(x,y,plane_array);
if get(handles.no_spb1_check, 'Value') == 1 || ...
        get(handles.no_plasmid_check, 'Value') == 1
    %set kmt1 to blank
    kmt1 = [];
    %display that kMT1 cannot be calculated
    set(handles.kmt1_text, 'String', 'NA');
else
    %calculate the distance in pixels from SPB1 to each plasmid pixel
    dist_array_spb1 = sqrt((plasmid_coords(:,1)*pixel_size - spb1(1)*pixel_size).^2 ...
        + (plasmid_coords(:,2)*pixel_size - spb1(2)*pixel_size).^2 ...
        + (plasmid_coords(:,3)*step_size - spb1(3)*step_size).^2);
    min_dist_spb1 = min(dist_array_spb1);
    kmt1 = min_dist_spb1;
    set(handles.kmt1_text, 'String', num2str(round(kmt1)));
end

if get(handles.no_spb2_check, 'Value') == 1 || ...
        get(handles.no_plasmid_check, 'Value') == 1
    %set kmt2 to blank
    kmt2 = [];
    %display that kMT1 cannot be calculated
    set(handles.kmt2_text, 'String', 'NA');
else
    %calculate the distance in pixel from SPB2 to each plasmid pixel
    dist_array_spb2 = sqrt((plasmid_coords(:,1) * pixel_size - spb2(1)*pixel_size).^2 ...
        + (plasmid_coords(:,2)*pixel_size - spb2(2)*pixel_size).^2 ...
        + (plasmid_coords(:,3)*step_size - spb2(3)*step_size).^2);
    min_dist_spb2 = min(dist_array_spb2);
    kmt2 = min_dist_spb2;
    set(handles.kmt2_text, 'String', num2str(round(kmt2)));
end


%if the plasmid is split then let the user choose which kmt to save
split_logical = get(handles.split_check,'Value');
if split_logical == 1
    kmt_array = get(handles.kmt1_check,'Value') -...
        get(handles.kmt2_check,'Value');
    switch kmt_array
        case 1
            kmt2 = [];
            set(handles.kmt2_text, 'String', 'NA')
        case -1
            kmt1 = [];
            set(handles.kmt1_text, 'String', 'NA')
    end
end

% add kmt 1 and 2 to handles structure for recording data
handles.kmt2 = kmt2;
handles.kmt1 = kmt1;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function calc_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function spb_1_button_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to spb_1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over spb_1_button.
function spb_1_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to spb_1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on spb_1_button and none of its controls.
function spb_1_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to spb_1_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function spb_2_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spb_2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function spb_2_button_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to spb_2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on spb_2_button and none of its controls.
function spb_2_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to spb_2_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function calc_button_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over calc_button.
function calc_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to calc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object deletion, before destroying properties.
function axes3_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function plane_slider_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to plane_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over plane_slider.
function plane_slider_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to plane_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on plane_slider and none of its controls.
function plane_slider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to plane_slider (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in record_button.
function record_button_Callback(hObject, eventdata, handles)
% hObject    handle to record_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This button will record the data in some sort of strucuture that
%can be saved later

%read in all the data we want
plasmid_plane = handles.plasmid_plane;
axis_length = handles.axis_length;
minor_axis = handles.minor_axis;
aspect_rat = handles.aspect_rat;
centroid = handles.centroid;
spb1 = handles.spb1;
spb2 = handles.spb2;
spindle_3D = handles.spindle_3D;
kmt1 = handles.kmt1;
kmt2 = handles.kmt2;
record_num = handles.record_num;
split_logical = get(handles.split_check,'Value');
im_check = handles.im_check;

%if the plasmid could not be dectected set kMTs and plasmid value to blank
if get(handles.no_plasmid_check, 'Value') == 1
    kmt1 = [];
    kmt2 = [];
    axis_length = [];
end

%if the plasmid is split, do not record spindle length
if get(handles.split_check, 'Value') == 1
    spindle_3D = [];
end

% put it all in the cell array that we created when we opened the image
handles.data_cell(record_num + 1,:) = {plasmid_plane, axis_length,...
    minor_axis, aspect_rat, centroid, spindle_3D, spb1, spb2,...
    kmt1, kmt2, split_logical, im_check};
%increase the record number
record_num = record_num + 1;
handles.record_num = record_num;
%update the record display text
set(handles.record_text,'String', num2str(handles.record_num));

%clear all the data
handles.plasmid_plane = [];
handles.axis_length = [];
handles.minor_axis = [];
handles.aspect_rat = [];
handles.spb1 = [];
handles.spb2 = [];
handles.spindle_3D = [];
handles.kmt1 = [];
handles.kmt2 = [];
handles.centroid = [];
handles.im_check = [];

%revert text boxes to original states
set(handles.axis_tag, 'String', 'Major Axis Length');
set(handles.minor_tag, 'String', 'Minor Axis Length');
set(handles.aspect_rat_tag, 'String', 'Aspect Ratio');
set(handles.kmt1_text, 'String', 'kMT 1');
set(handles.kmt2_text, 'String', 'kMT 2');
set(handles.spindle_nm, 'String', '3D Spindle');

%update the handles data structure
guidata(hObject,handles);


% --------------------------------------------------------------------
function save_data_Callback(hObject, eventdata, handles)
% hObject    handle to save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%load the cell array from the handles structure
data_cell = handles.data_cell;
uisave('data_cell');


% --- Executes on button press in no_spb1_check.
function no_spb1_check_Callback(hObject, eventdata, handles)
% hObject    handle to no_spb1_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_spb1_check


% --- Executes on button press in no_spb2_check.
function no_spb2_check_Callback(hObject, eventdata, handles)
% hObject    handle to no_spb2_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_spb2_check


% --- Executes on button press in clear_longaxis.
function clear_longaxis_Callback(hObject, eventdata, handles)
% hObject    handle to clear_longaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear the major axis length
handles.axis_length = [];
set(handles.axis_tag, 'String', 'Major Axis Length');
%clear the minor axis length
handles.minor_axis = [];
set(handles.minor_tag,'String', 'Minor Axis Length');
%clear the aspect ratio
handles.aspect_rat = [];
set(handles.aspect_rat_tag,'String','Aspect Ratio');
%clear the centroid measurement
handles.centroid = [];
guidata(hObject,handles);
%Change the figure to axes1, the GFP channel
subplot(handles.axes1);


% --- Executes on button press in clear_spindle.
function clear_spindle_Callback(hObject, eventdata, handles)
% hObject    handle to clear_spindle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Clear all spindle related info
handles.spb1 = [];
handles.spb2 = [];
handles.spindle_3D = [];
handles.kmt1 = [];
handles.kmt2 = [];

%reset all text boxes
set(handles.kmt1_text, 'String', 'kMT 1');
set(handles.kmt2_text, 'String', 'kMT 2');
set(handles.spindle_nm, 'String', '3D Spindle');


% --- Executes on button press in split_check.
function split_check_Callback(hObject, eventdata, handles)
% hObject    handle to split_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of split_check


% --- Executes on button press in no_plasmid_check.
function no_plasmid_check_Callback(hObject, eventdata, handles)
% hObject    handle to no_plasmid_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set the axis length in handles strucutre to blank


% --- Executes on button press in BG_button.
function BG_button_Callback(hObject, eventdata, handles)
% hObject    handle to BG_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the current plane number
slider_value = ceil(get(handles.plane_slider,'Value'));
%change the plot to axes1, the GFP image
subplot(handles.axes1);
%draw a rectangle
h = imrect;
%get position coordinates of rectangle
pos_bg = wait(h);
%sample the GFP image
im_mat_bg = handles.im_mat(round(pos_bg(2)):round(pos_bg(2))+round((pos_bg(4))),round(pos_bg(1)):...
    round(pos_bg(1))+ round(pos_bg(3)),slider_value);
%get average of that image
avg_bg = mean(im_mat_bg(:));
handles.avg_bg = avg_bg;
set(handles.BG_text, 'String', num2str(avg_bg));
guidata(hObject,handles);


% --- Executes on button press in imdilate_button.
function imdilate_button_Callback(hObject, eventdata, handles)
% hObject    handle to imdilate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%create a strucutral element
SE_disk = strel('disk', 1);

%load in the binary image
im_bin = handles.im_bin;

%apply the dilation
im_bin = imdilate(im_bin, SE_disk);

%show the new image
subplot(handles.axes2);
imshow(im_bin,[]);

%update the im_bin in the handles structure
handles.im_bin = im_bin;

%perform the same finishing operations as FG button
%load in the necessary variables from FG button
pos_fg = handles.pos_fg;

%pad the binary image back to the orginal size
[im_y, im_x, ~] = size(handles.im_mat);
x_pad_1 = zeros(pos_fg(4) + 1, pos_fg(1) - 1);
im_bin_resize = horzcat(x_pad_1,im_bin);
[~,im_bin_resize_width] = size(im_bin_resize);
x_pad_2 = zeros(pos_fg(4) +1, im_x - im_bin_resize_width);
im_bin_resize = horzcat(im_bin_resize, x_pad_2);

y_pad_1 = zeros(pos_fg(2) - 1, im_x);
im_bin_resize = vertcat(y_pad_1,im_bin_resize);
[im_bin_resize_height, ~] = size(im_bin_resize);
y_pad_2 = zeros(im_y - im_bin_resize_height, im_x);
im_bin_resize = vertcat(im_bin_resize, y_pad_2);
slider_value = ceil(get(handles.plane_slider,'Value'));
handles.plasmid_mask = im_bin_resize;
handles.plasmid_plane = slider_value;

%get the centroid of the resized binary image
if get(handles.no_plasmid_check,'Value') == 0
    centroid_struct = regionprops(im_bin_resize,'centroid');
    centroid = centroid_struct(1);
    handles.centroid = struct2array(centroid);
else
    handles.centroid = [];
    handles.plasmid_plane = [];
end
guidata(hObject,handles);


% --- Executes on button press in erode_button.
function erode_button_Callback(hObject, eventdata, handles)
% hObject    handle to erode_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%create a strucutral element
SE_disk = strel('disk', 1);

%load in the binary image
im_bin = handles.im_bin;

%apply the dilation
im_bin = imerode(im_bin, SE_disk);

%show the new image
subplot(handles.axes2);
imshow(im_bin,[]);

%update the im_bin in the handles structure
handles.im_bin = im_bin;

%perform the same finishing operations as FG button
%load in the necessary variables from FG button
pos_fg = handles.pos_fg;

%pad the binary image back to the orginal size
[im_y, im_x, ~] = size(handles.im_mat);
x_pad_1 = zeros(pos_fg(4) + 1, pos_fg(1) - 1);
im_bin_resize = horzcat(x_pad_1,im_bin);
[~,im_bin_resize_width] = size(im_bin_resize);
x_pad_2 = zeros(pos_fg(4) +1, im_x - im_bin_resize_width);
im_bin_resize = horzcat(im_bin_resize, x_pad_2);

y_pad_1 = zeros(pos_fg(2) - 1, im_x);
im_bin_resize = vertcat(y_pad_1,im_bin_resize);
[im_bin_resize_height, ~] = size(im_bin_resize);
y_pad_2 = zeros(im_y - im_bin_resize_height, im_x);
im_bin_resize = vertcat(im_bin_resize, y_pad_2);
slider_value = ceil(get(handles.plane_slider,'Value'));
handles.plasmid_mask = im_bin_resize;
handles.plasmid_plane = slider_value;

%get the centroid of the resized binary image
if get(handles.no_plasmid_check,'Value') == 0
    centroid_struct = regionprops(im_bin_resize,'centroid');
    centroid = centroid_struct(1);
    handles.centroid = struct2array(centroid);
else
    handles.centroid = [];
    handles.plasmid_plane = [];
end
guidata(hObject,handles);


% --- Executes on button press in kmt1_check.
function kmt1_check_Callback(hObject, eventdata, handles)
% hObject    handle to kmt1_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmt1_check


% --- Executes on button press in kmt2_check.
function kmt2_check_Callback(hObject, eventdata, handles)
% hObject    handle to kmt2_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmt2_check


% --- Executes on button press in clear_centroids_button.
function clear_centroids_button_Callback(hObject, eventdata, handles)
% hObject    handle to clear_centroids_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear the centroid_points
handles.centroid_points = [];

%clear the plot points by re-showing the im_mat in axes1 and axes2
slider_value = ceil(get(handles.plane_slider,'Value'));
subplot(handles.axes3);
imshow(handles.im_mat_2(:,:,slider_value),[]);
subplot(handles.axes1);
imshow(handles.im_mat(:,:,slider_value),[]);

%Update the handles structure
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
    case 'b' % select background
        BG_button_Callback(hObject, [], handles);
    case 'f' % select foreground
        FG_button_Callback(hObject, [], handles);
    case 'a' % get axis length
        major_axis_Callback(hObject, [], handles);
    case '1' % select first RFP
        spb_1_button_Callback(hObject, [], handles);
    case '2' % select second RFP
        spb_2_button_Callback(hObject, [], handles);
    case 'c'
        clear_centroids_button_Callback(hObject, [], handles);
    case 's' % calc spindles
        calc_button_Callback(hObject, [], handles);
    case 'r' % record data
        record_button_Callback(hObject, [], handles);
    case 'd' % dilate
        imdilate_button_Callback(hObject, [], handles);
    case 'e' % erode
        erode_button_Callback(hObject, [], handles);
end