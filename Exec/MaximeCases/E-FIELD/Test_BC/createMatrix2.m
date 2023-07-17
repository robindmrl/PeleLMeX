%name={'efieldx','efieldy','y_velocity','DriftFlux_H3O+_Y','temp'};
name={'x_velocity','y_velocity','temp'};

for i=1:length(name) %Loop for converting all files. Saves them automatically 
    % Erase all earlier variables from matlab.
    %clear all

    % Clear the plotting window.
    clc
    close all
    
    %Spliting the name to make it easier to change files name
    c1='*plt00039_'; 
    b=strcat(c1,name{1,i});
    c= split(b,'_');
    
    
    folder = '/home/maxdnz/PeleLMeX/Exec/MyCases/E-FIELD/Test_BC'; %input file folder
    fname = dir(fullfile(folder, strcat(b,'.mat')));  %Reads all .mat files

    % ref_ratios = input('Input ref ratios in [_,_,_ ...] for each level: ');
    ref_ratios = [2,2];
    ref=prod(ref_ratios);

    for f=1:numel(fname)
    % Read in the number of grids at each level, the locations and sizes 
    %  of those grids, and the data on each grid.
    [dim ngrid loc siz dat] = binread(fname(f).name);


    % Number of levels
    nlev = size(ngrid);

    % The size of the graph will be the size of the first level 0 grid.
    % If there is more than one level 0 grid, this will need to be modified.
    xlo = loc{1}{1}(1);
    xhi = loc{1}{1}(2);
    ylo = loc{1}{1}(3);
    yhi = loc{1}{1}(4);

    for n=1:ngrid(1)
      xlo = min(xlo,loc{1}{n}(1));
      xhi = max(xhi,loc{1}{n}(2));
      ylo = min(ylo,loc{1}{n}(3));
      yhi = max(yhi,loc{1}{n}(4));
    end

    %axis([xlo xhi ylo yhi]);

    ncont = 10;
    showBoxes = 1;


    %hold on;

     x_coarse=0;
     y_coarse=0;
    num_grids = ngrid(1); %First level
      for n=1:num_grids
          if loc{1}{n}(3)==0
            x_coarse =  x_coarse+size(dat{1}{n},2);
          end
          if loc{1}{n}(1)==0
             y_coarse = y_coarse+size(dat{1}{n},1);
          end
      end

    data(:,:,f) = zeros(y_coarse*ref,x_coarse*ref);  %creates big matrix to be filled

    for l=1:nlev
       num_grids = ngrid(l); 
       r=prod(ref_ratios(l:end));
       for n=1:num_grids
                data(round(loc{l}{n}(3)*y_coarse*ref/yhi+1):round(loc{l}{n}(4)*y_coarse*ref/yhi),round(loc{l}{n}(1)*x_coarse*ref/xhi+1):round(loc{l}{n}(2)*x_coarse*ref/xhi),f)=...
                 repelem(dat{l}{n},r,r);
       end
    end
    %   figure 
    %   h=pcolor(data(:,:,f))
    %   set(h, 'edgecolor','none')
    %   colormap(jet)
    %   pbaspect([xhi yhi 1])
    end
    
    save(strcat(c{2,1},'.mat'),'data');
    
    figure
      h=pcolor(data(:,:,f))
      set(h, 'edgecolor','none')
      colormap(jet)
      pbaspect([xhi yhi 1]) 
    fclose('all')   
end

