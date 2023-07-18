clear all

%% File info
main_folder = '/home/maxdnz/PeleLMeX/Exec/MyCases/E-FIELD/HPC_3.10/MATLAB_00000';
species = {'/Y(N2)','/Y(O2)','/Y(H2O)','/Y(OH)','Y(CO)','Y(CO2)','Y(CH4)','/Y(C)','/Y(CH)','/Y(CHOp)','/Y(H3Op)','/Y(O2n)','/Y(E)'};
Temp = {'/temp'};
Speed = {'/y'};
EF = {'/efieldx','/efieldy','/DriftFlux'};
MWval = {14.0067,31.999,18.01528,17.008,28.01,44.01,16.04,12.0107,13.01864,29.018,19.0232,31.999,0.00054858};

%% Import Data

% WARNING : ALL THE DATA IS IN THE WRONG WAY : Y 1ST THEN X

%{
% Species
for i=1:lenght(species)
    file = strcat(main_folder,species{1,i},'.mat');   

    if i==1
        A = zeros(lenght(species),size(importdata(file),1),size(importdata(file),2)); 
    end

    A(i,:,:) = importdata(file); %mass fractions
end
%}

% Temperature 
file2 = strcat(main_folder,Temp{1,1},'.mat');
B = zeros(size(importdata(file2),1),size(importdata(file2),2));
B(:,:) = importdata(file2);

% Velocity
file5 = strcat(main_folder,Speed{1,1},'.mat');
E = zeros(size(importdata(file5),1),size(importdata(file5),2));
E(:,:) = importdata(file5);

% Efield
for i=1:2
    file3 = strcat(main_folder,EF{1,i},'.mat');   

    if i==1
        C = zeros(2,size(importdata(file3),1),size(importdata(file3),2));
    end

    C(i,:,:) = importdata(file3);
end

%% Set up x,y axis
xmax = 18.5; %mm
ymax = 37; %mm
X = linspace(0,xmax,size(B,2));
Y = linspace(0,ymax,size(B,1));

%% Results check
%% Inflow profile - Speed
profile = E(1,:)*100;
plot(X,profile);

xlim([0,19])
title("Speed Low BC")
xlabel("x in mm")
ylabel("Speed in cm/s")

%% Low BC - Temp
TempBC = B(1,:);
plot(X,TempBC);

xlim([0,19])
title("Temp Low BC")
xlabel("x in mm")
ylabel("Temp in K")

%% Efield vectors
reducx = 10; %we keep 1 vector each reducx vectors along x axis
reducy = 40; %we keep 1 vector each reducy vectors along y axis

EFx = reshape(C(1,:,:),[size(C,2),size(C,3)]);
EFy = reshape(C(2,:,:),[size(C,2),size(C,3)]);
Norm_EF = sqrt(EFx.*EFx + EFy.*EFy);

X_red=zeros(1,floor(size(X,2)/reducx));
Y_red=zeros(1,floor(size(Y,2)/reducy));
EFx_reduce=zeros(floor(size(EFx,1)/reducy),floor(size(EFx,2)/reducx));
EFy_reduce=zeros(floor(size(EFy,1)/reducy),floor(size(EFy,2)/reducx));

for i = reducy:reducy:size(EFx,1)
    for j = reducx:reducx:size(EFx,2)
        EFx_reduce(floor(i/reducy),floor(j/reducx)) = EFx(i,j);
        EFy_reduce(floor(i/reducy),floor(j/reducx)) = EFy(i,j);
    end
end

for i = reducx:reducx:size(EFx,2)
   X_red(floor(i/reducx)) = X(i);
end

for j = reducy:reducy:size(EFx,1)
   Y_red(floor(j/reducy)) = Y(j);
end

quiver(X_red,Y_red,EFx_reduce,EFy_reduce,0.8,'r','linewidth',1.5); 

hold on
contour(X,Y,Norm_EF,'linewidth',1.5)
hold off

title("Efield vectors and EF norm")
xlabel("x in mm")
ylabel("y in mm")
legend('vectors','norm')

%% T profile in the flame - x = 30% of full size
height = 0.3;
ix = floor(size(B,1)*height);
TempFlame = B(ix,:);
plot(X,TempFlame);

xlim([0,19])
title("Temp profile at 30% height")
xlabel("x in mm")
ylabel("Temp in K")
