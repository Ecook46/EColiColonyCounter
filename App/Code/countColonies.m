function [RGB,objects] = countColonies(path)
%% Locate ecoli plate  
RGB = imread(path); %import image
grey = im2gray(RGB); % greyscale
bw = imbinarize(grey); % binarize

bw = imclearborder(bw); % remove regions touching the border

minsize = 5000; % set minsize for closing
bw = bwareaopen(bw, minsize,8); % closes whitespace holes from reflection

bw = bwmorph(bw, 'close', Inf); % performs an infinite number of dilation/erosions
bw = imfill(bw,'holes'); % fills in pixels not reached by background
se = strel('disk', 8); % create a structuring item
bw = imerode(bw, se); % erode image with the disk structures
bw = bwmorph(bw, 'majority', Inf); %sets pixel to 1 if 5 surrounding are 1
B = bwboundaries(bw,"noholes"); % get boundaries of objects in bw

stats = regionprops(bw,"Circularity","Centroid","Area"); % get circular properties 
threshold = 0.90; % threshold for circlularity
ratioArea = 0.70; % comparative area ratio 

remove = zeros(length(B)-1,1); %blank removal matrix
j = 1; % counter
mArea = max(cell2mat({stats(:).Area})); % find the largerst area
for k = 1:length(B)
  % index small and non circular regions for removal  
  if or(stats(k).Circularity < threshold,stats(k).Area < ratioArea * mArea) 
      remove(j) = k; % store index
      j = j+1; % increment
      continue
  end
end

for j = 1:length(remove) % remove invalid regions
    if remove(j) == 0 % prevent 0 indexing error
        remove(j) = [];
    end
end
 B(remove) = []; % remove invalid regions from boundaries
 stats(remove) = []; % remove invalid regions from stats

tempMax = 1; % temp max variable
% finding the rightmost valid region
for j = 1:length(stats)
    if stats(tempMax).Centroid(1) < stats(j).Centroid(1) % higher x = more right
        tempMax = j;
    end
end

  B = B(tempMax); % set B to rightmost valid region
  area = stats(tempMax).Area(1); % store area
  center = stats(tempMax).Centroid; % store centroid
 
  %% Isolate plate region

  % define boundaries of plate
  boundary = cell2mat(B(1));

  % set non plate region black
  [rows, columns] = size(bw); % get image size
  excludeMask = poly2mask(boundary(:,2),boundary(:,1),rows,columns); % create mask  
  RGB2 = imoverlay(RGB,~excludeMask,"Black"); % overlay the mask
  
  %% Isolate Colonies

  [bw,RGB2] = createMask(RGB2); % exclude objects not similar to ecoli

  % replacing black background with white
  whiteMask = repmat(all(~RGB2,3),[1 1 3]); 
  RGB2(whiteMask) = 255;
    
  %% Prepare for final clean and count

  % change to binary
  RGB2 = im2gray(RGB2);
  final = imbinarize(RGB2);

  % Crop image to smaller square
  radius = sqrt(area/pi)*1.3; % multiply by 1.3 for larger square bound
  % using trig to find corner points
  x1 = floor((radius*sind(45)) + center(1)); 
  y1 = floor((radius*cosd(45)) + center(2));
  x2 = floor((radius*sind(135)) + center(1));
  y2 = floor((radius*cosd(135)) + center(2));
  x3 = floor((radius*sind(225)) + center(1));

  final = final(y2:y1,x3:x2); % crop image by taking certain rows/columns
  RGB = RGB(y2:y1,x3:x2,:);

  final = ~final; % invert image

  
  % perform watershed transform to segment overlapping colonies
  dd = -bwdist(~final); % distance transform
  dd = imhmin(dd, 3); % avoid oversegmetation of shallow basins
  L = watershed(dd); % compute watersheds of the transform
  final(L == 0) = false; % apply the watershed

  se = strel('disk', 2); % create a structuring item
  final = imerode(final, se); % erode image with the disk structures

  final = bwareaopen(final,25); % close areas less than 25 pixels

  %% Count final objects

  objects = regionprops(final,"Centroid"); % get all remaining objects centroids
  
end
