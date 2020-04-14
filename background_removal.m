%                    Background Removal Function                          %
%-------------------------------------------------------------------------%
% Goal: Remove background images from holograms
%-------------------------------------------------------------------------%
% Use: Place monochrome background images in bg*.png format into the
%           ./background_images/ folder, with * corresponding to image
%           number.
%      Place hologram to be processed in the working directory of the
%           script.
%-------------------------------------------------------------------------%
% Outputs: The processed hologram will be saved as a .mat file in the
%          working directory; Function will pass background-removed image.
%-------------------------------------------------------------------------%

function im_divided = background_removal(input_image)

% Set hologram file name here:

    holo_image = double(imread(input_image));
    
figure(1); imshow(holo_image, []); title('Original Image')


    
% Get list of background files

    bgfiles=dir(fullfile('./background_images/','bg*.png'));
    num_bgs=length(bgfiles);

    
    
% Initialize bg_sum variable

    bg_sum = 0;

    
    
% Compute & display average of background

    for i=1:num_bgs
        a=double(imread(strcat('./background_images/',bgfiles(i).name)));
        bg_sum=bg_sum+a;
    end

    bg_avg = bg_sum / num_bgs;

    figure(2); imshow(bg_avg,[]); title('Averaged Background')
    


% Compute & display processed hologram

    % Should I divide by the intensity, or the square-root of intensity
    % (the amplitude) here? dividing by intensity looks better, but their
    % formula uses square-root. Need to research this.
    im_divided = (holo_image-bg_avg)./bg_avg;

    figure(3); imshow(im_divided,[]); title('Background Removed')


    
% Save hologram

    save('holo_processed.mat','im_divided');

