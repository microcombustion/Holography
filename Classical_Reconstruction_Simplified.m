%% Load Hologram

% hol_input = double(imread('24mm_ruler.png'));
% load('holo_processed.mat');


% hol_input = im_divided;
% figure(1); imshow(hol_input,[]); title('Input Hologram')

hol_input = background_removal('24mm_ruler.png');

%% Reconstruction Z-Range

    lb      = 0.069   ;	% Lower Bound
    stepZ   = 0.0005  ;	% Step
    ub      = 0.071   ;	% Upper Bound

    Z=lb:stepZ(1):ub  ; % Z-Range
%    Z(Z==0) = eps     ; % eps is spacing of floating point numbers, 2^(-52)
    


%% Camera & Laser Information

    K       = 0.8367  ; % Fill Factor
    pix     = 3.6e-6  ;	% Pixel Size
    lambda  = 532e-9  ;	% Laser Wavelength



%% WHY THIS?    
% Mean intensity for all pixels for the input hologram is subtracted
% from the input.

%     hol_input = hol_input-mean(hol_input(:));


%% Zero Pad Hologram for FFT

    hol_size = size(hol_input);

    holo_pad = zeros(hol_size*2);
    holo_pad( hol_size(1)/2+1 : 3*hol_size(1)/2 , hol_size(2)/2+1 : 3*hol_size(2)/2) = hol_input;
    
    figure(4); imshow(holo_pad,[]); title('Input Hologram for Reconstruction, with Zero-Padding')



%% Fresnel Transform

Reconst     = zeros(hol_size(1),hol_size(2),numel(Z));  % Initialize
hol_sizep   = size(holo_pad);                           % Size of Padded Hologram
dum         = zeros(hol_sizep); 

% Nyquist frequency
    % Calculated from shorter side. 
    % Something like (1024*pix^2 / (lambda)) * 2
    % Need to understand this.
Z_nyq = min(hol_sizep*pix)^2 / (min(hol_sizep)*lambda); 

c = 0;  % Counter
    
% Form meshgrid 
    % p is x-range -1280 to 1279 (2560 elements)
    % j is y-range -1024 to 1023 (2048 elements)
[p, j] = meshgrid(floor(-hol_sizep(2)/2):floor(hol_sizep(2)/2)-1,...
                  floor(-hol_sizep(1)/2):floor(hol_sizep(1)/2)-1);

              
    
%%% Reconstruction using Transfer function frequency %%%

% Looking for Z positions <= to Nyquist frequency
    % This is used to determine reconstruction method.
    % For Z positions <= Z_nyq, reconstruction uses transfer
    % function frequency. For positions > Z_nyq, reconstruction uses
    % impulse response frequency.
    ind = find(abs(Z) <= Z_nyq);

if numel(ind)>0
    Z_bi = Z(ind);

    % Polar Conversion. Zero is at [1025 1281].
        % Need to work out how this differs from impulse-response radius
    Ro = (j/hol_sizep(1)/pix).^2+(p/hol_sizep(2)/pix).^2;
      
    for ii = 1:numel(Z_bi); c=c+1;

        Hz =  exp(1i*pi*lambda*Z_bi(ii)*Ro);

        dum = ifft2(fft2(holo_pad).*fftshift(Hz)); % Convolution

        Reconst(:,:,c) = dum(floor(hol_size(1)/2)+1:3*floor(hol_size(1)/2),...
                             floor(hol_size(2)/2)+1:3*floor(hol_size(2)/2));

        % What's this for? Why did they comment it out?
        % Reconst(:,:,c)=(Reconst(:,:,c)-min(min(Reconst(:,:,c))))./(max(max(Reconst(:,:,c)))-min(min(Reconst(:,:,c))));

        fprintf('Reconstruction using Transfer Function Hz for z=%g is done\n',Z_bi(ii));
    end
end    



%%% Reconstruction using impulse-response frequency %%%
ind = find(abs(Z)>Z_nyq); % in space domain

if numel(ind)>0
   Z_sm = Z(ind);
   
   Ro = (j*pix).^2+(p*pix).^2;

    for ii = 1:numel(Z_sm); c=c+1;
        
        hz = 1/1i/lambda/Z_sm(ii) .* exp(pi*1i*Ro/lambda/Z_sm(ii));

        dum = ifftshift( ifft2( fft2(holo_pad) .* fft2(hz)) ); % Convolution

        Reconst(:,:,c) = dum(floor(hol_size(1)/2)+1:3*floor(hol_size(1)/2),...
                             floor(hol_size(2)/2)+1:3*floor(hol_size(2)/2));

        % Reconst(:,:,c)=(Reconst(:,:,c)-min(min(Reconst(:,:,c))))./(max(max(Reconst(:,:,c)))-min(min(Reconst(:,:,c))));

        fprintf('Reconstruction using Impulse Response hz for z=%g is done\n',Z_sm(ii));

    end
end

titlestring = sprintf('Reconstruction @ Z = %3.3f m',lb);

figure(5); imshow(real(Reconst(:,:,1)),[]); title(titlestring)
recon = abs(Reconst(:,:,1));
figure(6); imshow(recon,[]); title(titlestring);