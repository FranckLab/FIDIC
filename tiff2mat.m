function [cellIMG,filename] = tiff2mat(Folder,ext,s)
    %Read images and write them out in .mat

    %Load all of the files
    files = dir(strcat(Folder,'\*',ext));

    %Determine the number of files
    if nargin<3
        s = length(files);
    end
    
    % Loop to open files- reads them in in alpha-numeric order
    for ii = 1:s
        READ = imread(strcat(Folder,'\',files(ii).name));
        IMG(:,:,ii) = double(READ(:,:,1));

    % Option to plot the images

%         imshow(IMG(:,:,ii))
%         drawnow
    end

    cellIMG = cell(1);

    for ii = 1:s
        cellIMG{1} = IMG(:,:,ii);
        filename = strcat('IDIC_image',num2str(ii-1));
        save(filename,'cellIMG');
    end
    
    filename = 'IDIC_image*';
    
end