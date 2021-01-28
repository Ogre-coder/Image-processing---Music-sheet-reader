%1. Read a music sheet 

music_sheet= im2double(rgb2gray(imread('Butyandthebeast.png')));

%2. break the sheet into sections based on staves 

%a)Increase the contrast of music sheet to be black and white only 

Threshold = graythresh(music_sheet);
Black_White = imbinarize(music_sheet,Threshold);

figure
imshow(Black_White)

%b) Make the black note appear to be white note on a black ground
%   so lines can be read using erosion

invert_music_sheet = ones(size(music_sheet)) - Black_White;
music_sheet = invert_music_sheet;
Original_Page = ones(size(music_sheet)) - music_sheet;

%c) Erode the staff lines uisng window of 65 pixels to filter out
%   other elements and keep only staves 

staves_only = imerode(music_sheet, ones(1,65));  
music_sheet = staves_only;
figure
imshow(music_sheet)

%d) Enhance the white color of stafflines using dilation
dilated_image = imdilate(music_sheet, ones(1,65));
music_sheet = dilated_image;
figure
imshow(music_sheet)

%e) Combine the length of entire staff lines into 1 column in order to 
%   identify which rows are staff lines 

staff = sum(music_sheet,2)/size(music_sheet,2);

%f) Increase contrast of black and white on the horizontal sum of the image
%   to get rows of staff lines
lines = im2bw(staff,graythresh(staff));

% g) Find an empty point to break the image into segments

%    i)Make area of breakpoints known as black 
     blank_area = zeros(size(lines));
     
%    ii)Switch the white staff lines to black lines  
     lines_to_black = ones(size(lines)) - lines;
     
%    iii)Scanning the music sheet based on the 8 pixel gaps between each staff line
     for i = 8:round(size(music_sheet,1)/8) 
         
%       iv)make a structuring element to color black lines white, as well 
%          as the gap between when found      
        SE = ones(i,8);
        
%       v)mark staff lines area semi-white or white and combine them with black space; 
%       black spaces are now parts where we can break the music sheet 
        blank_area = blank_area + imerode(lines_to_black, SE);
     end
    
% h) threshold set up using the median filtering snippet found online 
%    from stack exchange
median_filter = blank_area .* (blank_area > median(blank_area));
    
%    i) discard zeros and take 10th percentile    
if(sum(median_filter) > 0)
    median_filter = sort(median_filter((median_filter > 0)));        
    median_filter = median_filter(round(.1* length(median_filter)));
else
    % ii) remove the entries with the maximum to get the threshold level
    divided = blank_area;
    divided(divided > 0.9*max(blank_area)) = [];
    median_filter = median(divided);
end

% i) Perform thresholding to increase the contrast of v) lines detection 
blank_area = blank_area > median_filter;

% j) use bwconncomp to find the number of segments and the size of cell that needed to break 
%    from the markings of v) we have marked . 
Per_segment = bwconncomp(1-blank_area);

% k) make a divider template 
divider_template = zeros(size(staff));

% l) make a cell to store the staves
staves = cell(Per_segment.NumObjects,1);    
    for i= 1:Per_segment.NumObjects
        each_cell = Per_segment.PixelIdxList{i};
        %{
        Threshold = round((max(Per_segment.PixelIdxList{i}) - min(Per_segment.PixelIdxList{i}))/4);
        extended_cells = min(Per_segment.PixelIdxList{i}) - Threshold : max(Per_segment.PixelIdxList{i}) + Threshold;
       
        divider_template(extended_cells(1)) = 1;
        divider_template(extended_cells(end)) = 1;
        %}
        staves{i} = Original_Page(Per_segment.PixelIdxList{i},:);
    end  
% m) output the staves
figure 
imshow(staves{1,1})
figure
imshow(staves{2,1})
figure
imshow(staves{3,1})
figure
imshow(staves{4,1})
        
% 3. Get the row numbers for the staff lines

%   a) make template for horizontal length of staff line 
%      and vertical length of 5 staff lines
    horizontal_length = zeros(length(staves),1);
    max_length = zeros(length(staves),1);
    staff = zeros(5,length(staves));
%   b) scan thru each section of 5 staves lines
    for i=1:length(staves)
%       i) sum vertically (only between the top and bottom staff line) to get the start of the line
        sum_5lines = sum(ones(size(staves{i})) - staves{i},2);
        
%    c) scan thru the gap and midpoint until we have 5 lines that match the
%           stave
        for gap = 8:length(sum_5lines)/5
            for midpoint = 2 : length(sum_5lines)/2;
 %              i) case1 scans thru the region within the line
                case2 = (midpoint:gap:(4*gap + midpoint));
 %              ii)  case2 scans thru the region above the line
                case2 = (midpoint+1:gap:(4*gap + midpoint+1));
 %              iii)    case3 scans thry the region below the line              
                case3 = (midpoint-1:gap:(4*gap + midpoint-1));
 %   d) if not reached the end of the staff lines, sum up all the cases to
 %      take account of the little area above and below the line as a
 %      safety cushion to capture the whole image 
                if (case2(5) < length(sum_5lines))
                    horizontal_length(i) = sum(sum_5lines([case2 case2 case3]));
                end
 %    e)if reached the end of the line, store the row number of staff lines 
 %       into matrix staff
                if (horizontal_length(i) > max_length(i))
                    max_length(i) = horizontal_length(i);
                    % store the result
                    staff(:,i) = case2;
                end
            end
        end
    end
    
% f) Output the row location of each staff line    
staff

%4) Find the key signature that signals the start of the music sheet

%a) Make a cell that contains the keys 

keys = cell(size(staves));

%b) Make a divider to segment the keys 

dividers = zeros(size(staves));

%c) Find the start of the music sheet 
    for i=1:length(staves)
        
%       i) Erase the bar lines    
        case2 = staff(:,i) - 1; % case2: 1 below the 5 lines 
        case3 = staff(:,i) + 1; % case3: 1 above the 5 lines 
        
%       ii) Output staves that without bars        
        barless = staves{i}; 
        barless([case2 staff(:,i) case3]) = 1;
        spacing = round(mean(mean(diff(staff,1),1)));
        sum_5lines = sum(~barless(staff(1,i): staff(5,i), :), 1);
        
        T = quantile(sum_5lines, .3);
        start0 = find(sum_5lines > T, 1);
        
        clef = ~staves{i}(staff(1,i):staff(5,i), start0:start0 + 3*spacing);
        com = cumsum(sum(clef,1) / sum(sum(clef)));
        start1 = find(com > 0.5, 1);
        keys{i} = ~staves{i}(:, start0+start1:end);
        
        new_start = start0 + start1;
            % check for activity on both size
        for j=3:size(keys{i},2)-3
            for k=1:5
                window = [0 3 5 3 0;...
                          0 2 4 2 0;...
                          0 1 3 1 0;...
                          0 0 2 0 0;...
                          0 0 1 0 0];  

                activity_above = sum(sum(window .* keys{i}(staff(k,i)-4:staff(k,i),j-2:j+2)));
                activity_below = sum(sum(flipud(window) .* keys{i}(staff(k,i):staff(k,i)+4,j-2:j+2)));
                T = 15;
                % Erase it if there isn't enough activity
                if ((activity_above < T) && (activity_below < T))
                    barless(staff(k,i)-4:staff(k,i)+4,j) = zeros(9,1);
                end
            end
        end
        sum_5lines = sum(barless,1);
%         (staff(1,i): staff(5,i))
%         collapsed = collapsed(new_start:round(.75* size(staves{i},2)));
        
%         T = quantile(collapsed, 0.01);
        start2 = find(sum_5lines == 0);
%         start2 = find(collapsed < T);
        if (length(start2) < 1)
            start2 = round(spacing/2);
        else
            start2 = start2(1);
        end
        dividers(i) = start0 + start1 + start2;
        keys{i} = keys{i}(:, start2:start2 + 10*spacing);
         imshow(~keys{i});
%         pause;        

    end
    
    figure
    imshow(~keys{1,1})
    figure
    imshow(~keys{2,1})
    figure
    imshow(~keys{3,1})
    figure
    imshow(~keys{4,1})
    
    size(staves)
    staves{1}
    solid_notes = cell(size(staves));
    for i=1:length(staves)
        spacing = round(mean(mean(diff(staff,1),1)));
        solid_notes{i} = ~staves{i};
    
     % erode the image using a small circle
        SE = strel('disk', round(spacing/4));
        solid_notes{i} = imerode(solid_notes{i}, SE);
        solid_notes{i} = imdilate(solid_notes{i}, SE);
       
        % Look at connected components
        CC = bwconncomp(solid_notes{i});
        stats = regionprops(CC, 'Eccentricity', 'Area', 'Centroid');
        objs = [];
        areas = [];
        matches = zeros(CC.ImageSize);
        
        for j=1:CC.NumObjects
            % Reject components that are less than half the note area
%             if (stats(j).Area < 0.125*spacing*spacing)
%                 continue;
%             end
            
            % Keep the slghtly eccentric circles that aren't too close to
%            the clef
            if (stats(j).Eccentricity < 0.8 && ...
                    stats(j).Eccentricity > .6 && ...
                    stats(j).Centroid(1) > dividers(i)) % has to bee past the clef
                objs = cat(1, objs, j);
                areas = [areas ; stats(j).Area];
            end
        end
        
         if (std(areas) > 0.25*spacing)
            objs = objs(areas > (median(areas) - 1.75*std(areas)));
        end
        
        solid_notes{i} = zeros(length(objs),2);
        for j=1:length(objs)
            matches(CC.PixelIdxList{objs(j)}) = 1;
            solid_notes{i}(j,:) = stats(objs(j)).Centroid;
        end
        % Flip so that the first item is the y coordinate (row)
        % and the second is the x coordinate (col)
        solid_notes{i} = fliplr(solid_notes{i});
        
        r = double(staves{i} | matches);
        g = double(staves{i} & ~matches);
        b = double(staves{i} & ~matches);
        
        colorimg = cat(3, r,g,b);
        figure
         imshow(colorimg);
    end