clear all

%% File info
main_folder = '/home/maxdnz/PeleLMeX/Exec/MyCases/E-FIELD/Test_BC';
Temp = {'/temp'};
Speed = {'/x','/y'};

%% Import Data

% WARNING : ALL THE DATA IS IN THE WRONG WAY : Y 1ST THEN X

% Temperature 
file = strcat(main_folder,Temp{1,1},'.mat');
A = zeros(size(importdata(file),1),size(importdata(file),2));
A(:,:) = importdata(file);

% Velocity
for i=1:2
    file2 = strcat(main_folder,Speed{1,i},'.mat');   

    if i==1
        B = zeros(2,size(importdata(file2),1),size(importdata(file2),2));
    end

    B(i,:,:) = importdata(file2);
end

Vx = reshape(B(1,:,:),[size(B,2),size(B,3)]);
Vy = reshape(B(2,:,:),[size(B,2),size(B,3)]);
V_norm = sqrt(Vx.*Vx + Vy.*Vy);

%% Set up x,y axis
xmax = 18.5; %mm
ymax = 37; %mm
X = linspace(0,xmax,size(A,2));
Y = linspace(0,ymax,size(A,1));

%% Initial condtions check
%% Inflow profile - Speed
profile = Vy(1,:)*100;
plot(X,profile);

xlim([0,19])
title("Speed Low BC")
xlabel("x in mm")
ylabel("Speed in cm/s")

%% Low BC - Temp
TempBC = A(1,:);
plot(X,TempBC);

xlim([0,19])
title("Temp Low BC")
xlabel("x in mm")
ylabel("Temp in K")

%% Blob Temp profile - y = 0.002, ymax = 0.037
ymax = 0.037;
yb = 0.002;
ib = floor(size(A,1)*yb/ymax);

BlobBC = A(ib,:);
plot(X,BlobBC);

xlim([0,2.5])
title("Blob Temp profile")
xlabel("x in mm")
ylabel("Temp in K")
