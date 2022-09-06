%% GENERATE SYNTHETIC MICROSTRUCTURE
% This code was modified so that the user can hard input the aspect ratio
% values, a and b, to keep the particle size consistent with volume
% fraction.

% First, clear all variables.

% Initialize variables
Vf = 0; 
nparticles = 0;
reject = 0;
noverlap = 0;

Vf_list=[];
for i=1:1
%% Input Microstructure Parameters
    Vf_max = 0.01;  %+(0.5-0.02)*rand;
    Vf_list(i) = Vf_max;
    a0 = 40;
    b0 = 40;
    image_size = 1024;

    if a0 < b0, c0 = a0; a0 = b0; b0 = c0; clear c0; end
    I_ellipse = image_ellipse(a0, b0, 0);


        N = ceil(image_size^2*Vf_max/sum(I_ellipse(:)));
        ratio = a0/b0; % a/b ratio
        diam = b0;



%% Select Particle Size Distribution
    N = round(1.15 * N); % Just in case, 15% cushion in particles
    List(1:N,1) = a0; List(1:N,2) = b0;
        

%% Select Particle Orientation Distribution



    Theta_List = linspace(0,pi,N+1);
    Theta_List(N+1) = [];


%% Determine Nearest Neighbor Parameters
% For ease of finding nearest neighbors for ellipse overlap criterion.

dltmp = diam;
nglnk = ceil( image_size / dltmp );
dlnk = image_size / nglnk;
lnk = zeros(nglnk,nglnk);
ndxy = 5 + ceil(2*a0/dlnk);

%% Input Filename

% Save image as filename
string = 'fiber_vf';
string_path = pwd;
string = strrep(string,'_vf',sprintf('_%.4f',Vf_max)); 
filename = sprintf('%s.tif',string);



%% Microstructure Generation Routine

Rej(1:N) = 0;
circles = zeros(N,5);
Image = zeros(image_size,image_size);
Image_elements = numel(Image);
Image_sum = 0;

tic
while Image_sum/Image_elements < Vf_max

    % First, randomly select the ellipse center coordinates.  Assign the
    % 'a' and 'b' parameter based on the largest particles first.  The
    % 'theta' parameter is randomly selected.

    x0 = ceil(rand*image_size); 
    y0 = ceil(rand*image_size);
    if (nparticles+1) <= size(List,1)
        a0 = List(nparticles+1,1);
        b0 = List(nparticles+1,2);
    end
    theta_rand = ceil(rand*length(Theta_List));
    theta0 = Theta_List(theta_rand);

    % If this is the first particle, no need to check for overlapping
    % ellipses.  Otherwise, check this potential ellipse against all other
    % ellipses in the neighborhood and return '0' if none overlap & '1' if
    % one overlaps this ellipse.

    if nparticles > 0;
        [overlap] = check_overlap_ellipses(lnk,circles,dlnk,nglnk,ndxy,...
            image_size,x0,y0,a0,b0,theta0);
    end

    % If there is no overlap, then keep the ellipse; otherwise, discard it.

    if nparticles == 0 || overlap == 0;

        % First, advance the number of particles and store all relevant
        % information in the matrix 'circles'.  Then, store the ellipse
        % center in the matrix 'lnk' for easy finding of neighboring
        % ellipses for overlap checks.  Finally, store the rejection
        % numbers and remove the theta from the theta dsitribution list.

        nparticles = nparticles + 1;
        circles(nparticles,1:5) = [x0 y0 a0 b0 theta0];
        [lnk] = store_circle(lnk,x0,y0,dlnk,nparticles);
        Rej(nparticles) = reject;
        reject = 0;
        Theta_List(theta_rand) = [];

        % Create the pixelated ellipse in matrix 'I_ellipse'.  This also
        % entails determining the size of 'I_ellipse' and the low and high 
        % bounds for inserting this ellipse into the main image matrix, 'Image'.

        if nparticles == 1
            I_ellipse = image_ellipse(a0, b0, theta0);
        end
        [nx,ny] = size(I_ellipse);
        if mod(nx,2)==1; xlo = (nx + 1)/2-1; xhi = xlo; end
        if mod(ny,2)==1; ylo = (ny + 1)/2-1; yhi = ylo; end
        if mod(nx,2)==0; xlo = (nx)/2; xhi = xlo-1; end
        if mod(ny,2)==0; ylo = (ny)/2; yhi = ylo-1; end

        % Calculate the correct bounds for inserting 'I_ellipse'.  This
        % includes wrapping the bounds through the periodic boundaries.

        nlo = x0 - xlo; nhi = x0 + xhi;
        ix = mod(nlo:nhi,length(Image));ixt = ix==0; ix(ixt)=length(Image);
        nlo = y0 - ylo; nhi = y0 + yhi;
        iy = mod(nlo:nhi,length(Image));iyt = iy==0; iy(iyt)=length(Image);

        % Now, add 'I_ellipse' to 'Image'.  If this ellipse overlaps with
        % another ellipse, then the binary matrix 'Image' will have a
        % maximum value greater than 1.  If this is the case, then
        % 'I_ellipse' will be moved until it is not overlapping.  This
        % seldom occurs, though, and is mainly dealt with in the symbolic
        % step of this algorithm.

        Image(ix,iy) = Image(ix,iy) + I_ellipse;
        if max(max(Image(ix,iy)))>1
            noverlap = noverlap + 1;
            [Image,x0,y0] = move_circle2(Image,I_ellipse,xlo,xhi,ylo,yhi,ix,iy);
            circles(nparticles,1) = x0; circles(nparticles,2) = y0;
            a = lnk == nparticles; lnk(a) = 0;
            [lnk] = store_circle(lnk,x0,y0,dlnk,nparticles);
        end

        Image_sum = Image_sum + sum(I_ellipse(:));

    else

        % If the ellipse overlaps, then add 1 to the reject counter and
        % start the whole process over again.

        reject = reject + 1;

    end
end
disp(['Symbolic generation of synthetic microstructure (sec): ' num2str(toc)])

%% Create/Save the image

tic
%figure;
%imshow(Image);


%imwrite(Image,filename,'tif'); 
%pathname = fileparts('E:\OneDrive - Northwestern University\Liu Research\SCA related\mesh generation\RVE_2D\RVE_gen');
Vf_actual=sum(Image(:))/numel(Image);
name = sprintf('Vf_%.4f.mat',Vf_max);
%name = fullfile(pathname,name);
save(name,'Image');
%disp(['Image viewing and saving (sec): ' num2str(toc)]);

%% Initial Ellipse Parameters Comparison

% This section compares the initial area fraction, number of particles, and
% a/b ratio with these metrics after pixelating the ellipses.

%disp('*** Computed Image Microstructure Metrics ***');
%disp(['GOAL - Area fraction: ' num2str(Vf_max)]);
disp(['Actual pixelated area fraction: ' num2str(Vf_actual)]);
%disp(['GOAL - Number of particles: ' num2str(N)]);
%disp(['Actual number of particles: ' num2str(nparticles)]);

%% Save Data to Excel 
% Save the ellipse parameters to Excel for future retrieval and
% postprocessing, including the nearest neighbor distribution.


    
%    a = circles(:,3)~=0; circles = circles(a,:);    
%    Excel_fileName = sprintf('%s\\ellipses.xls',string_path);
%    ColHeaders = {'x0','y0','a0','b0','theta0'};
%    sheet = strrep(string,'_vf',sprintf('_%dvf',Vf_max*100));
%    dist = 1:length(circles(:,1));
%    xlswrite(Excel_fileName, ColHeaders, sheet, 'A1');
%    xlswrite(Excel_fileName, circles, sheet, 'A2');
%    xlswrite(Excel_fileName, circles(:,5)/pi*180, sheet, 'E2');
%    deleteEmptyExcelSheets(Excel_fileName);


% End of the iteration loop for Vf_max
end
fileID = fopen('id_set.txt','w');
fprintf(fileID,'%.4f\n',Vf_list);
fclose(fileID);