close all;
clear all;
clc;
% Load the image
img = imread('image.jpg');  % Use the path to your image
gray_img = rgb2gray(img);   % Convert to grayscale

% Display the image
imshow(gray_img);
title('Select the first lane manually');

% Manually draw the first rectangle
first_rect = drawrectangle();  
rect_position = round(first_rect.Position);  % Get position of the first lane [x, y, width, height]

% Set the width and height based on the first lane
lane_width = rect_position(3);
lane_height = rect_position(4);

% Initialize an array to store positions of the ROIs
rect_positions = rect_position;  % Start with the first ROI

% Initialize a cell array to store the signal profiles
signal_profiles = {};

% Loop to allow the user to click the top-left corner for more rectangles
continue_drawing = true;
lane_index = 1;

while continue_drawing
    % Ask if the user wants to add another rectangle
    answer = questdlg('Do you want to add another rectangle by clicking the top-left corner?', ...
                      'Continue Drawing', 'Yes', 'No', 'Yes');
    
    if strcmp(answer, 'Yes')
        % Display the image and let the user click to define the top-left corner
        imshow(gray_img);
        hold on;  % Keep the image displayed for drawing
        title('Click the top-left point for the next rectangle');
        
        % Let the user click on the image for the top-left corner
        [x_coord, y_coord] = ginput(1);
        
        % Define the new rectangle position based on the clicked point and same size as the first rectangle
        new_rect_position = [round(x_coord), round(y_coord), lane_width, lane_height];
        
        % Store the new rectangle position
        rect_positions = [rect_positions; new_rect_position];
        
        % Draw the new rectangle in real-time on the image
        rectangle('Position', new_rect_position, 'EdgeColor', 'r', 'LineWidth', 2);
        hold off;  % Update the image to show the newly drawn rectangle
        
        % Crop the ROI for the current lane
        roi = imcrop(gray_img, new_rect_position);
        
        % Sum intensity vertically (along columns)
        signal_profile = sum(roi, 2);
        
        % Store the signal profile in the cell array
        signal_profiles{lane_index} = signal_profile;
        lane_index = lane_index + 1;
        
    else
        % Stop the loop if the user chooses "No"
        continue_drawing = false;
    end
end

% Save signal profiles to Excel
num_profiles = length(signal_profiles);
max_length = max(cellfun(@length, signal_profiles));  % Find the longest profile

% Initialize a table to store the data
signal_table = table;

% Loop through all profiles and store them in the table
for i = 1:num_profiles
    profile_length = length(signal_profiles{i});
    % Pad shorter profiles with NaN to match the length of the longest profile
    padded_profile = [signal_profiles{i}; nan(max_length - profile_length, 1)];
    signal_table = [signal_table, array2table(padded_profile, 'VariableNames', {sprintf('Lane_%d', i)})];
end

% Write the table to an Excel file
writetable(signal_table, 'signal_profiles.xlsx');

% Final display of all selected lanes on the image
imshow(gray_img);
title('Selected ROIs');
hold on;
for i = 1:size(rect_positions, 1)
    rectangle('Position', rect_positions(i,:), 'EdgeColor', 'r', 'LineWidth', 2);  % Display all lanes
end
hold off;

% Initialize the figure to show profiles
figure;
hold on;
for i = 1:num_profiles
    % Plot the signal
    plot(signal_profiles{i}, 'DisplayName', sprintf('Lane %d', i));
end

title('Signal Intensity Profiles');
xlabel('Position (pixels)');
ylabel('Intensity');
legend show;
