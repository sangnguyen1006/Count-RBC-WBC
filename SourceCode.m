%% RBC, WBC
%% Read image, convert RGB -->GRAY,HSV
[filename, pathname]=uigetfile({'*.jpg;*.png;*.dcm;*.bmp;*.gif'});
RGB=imread(strcat(pathname,filename));

if size(RGB,3)==3
    GRAY= rgb2gray(RGB);
    HSV= rgb2hsv(RGB);
end

%% Convert to binarize image
GRAY=adapthisteq(GRAY);
S=HSV(:,:,2);                             % S-Channel 

GRAYb= imbinarize(GRAY,graythresh(GRAY)); % Global image threshold using Otsu's method
Sb = imbinarize(S,graythresh(S));         % Global image threshold using Otsu's method
se = strel('disk',30);                    % Specify a radius of 5 pixels so that the largest gap gets filled.
Sb= imclose(Sb,se);  

%% Segmentation WBC,Plateles
Platelets= bwareafilt(Sb,[0 1500]);       % take area (0-1000 pixels)
WBC =bwareaopen(Sb,1500);                 % removal Platelets (area<1000 pixels)

%% Segmentation RBC
%preRBC=abs(Sb-(~GRAYb)); % removal WBCs
preRBC=~GRAYb-Sb;
for i=1:size(preRBC,1)
    for j=1:size(preRBC,2)
        if preRBC(i,j)<0
            preRBC(i,j)=0;
        end
    end
end

RBC= bwareaopen(preRBC,400);              % reducing noise(area<400 pixels)
se = strel('disk',5);                     % Specify a radius of 5 pixels so that the largest gap gets filled.
RBC= imclose(RBC,se);                     % smooth with the shape of disk

RBC= imfill(RBC,'holes');                 % filling holses
%RBC =bwareaopen(RBC,1000);               % removal Platelets (area<1000 pixels)

%% Hough transform method
% Rmax, Rmin 
[B,L] = bwboundaries(RBC,'noholes');      % tim duong bien
STATS = regionprops(L, 'all');
C0=zeros(size(B,1),1);                    % khoi tao ma tran 0
for i=1:size(B,1) 
    C0(i)=(STATS(i).Perimeter)^2/STATS(i).Area; % do nen
end
number=1:1:length(C0);
for i=1:length(C0)-1                      % sort
    for j=i+1:length(C0)
        if C0(j)<C0(i)
            tam=C0(j);
            C0(j)=C0(i);
            C0(i)=tam; 
            
            tam=number(j);
            number(j)=number(i);
            number(i)=tam;
        end
    end
end
for i=1:length(C0)                        % xac dinh hinh dang tron 14<=C0<15
    if C0(i)>=15
        mark=i-1;
        break;
    end
end
R=zeros(mark,1);
for i=1:mark
    R(i)=(2*STATS(number(i)).Area)/STATS(number(i)).Perimeter; % tinh R
end
R=sort(R);

Rmax=R(length(R));
Raverage=sum(R)/length(R);
Rmin=round(Raverage-Raverage*0.3);
Rmax=round(Rmax+Rmax*0.2+(Rmax-Raverage)*0.2);
%Hough transform 
[center, radius, metric] = imfindcircles(RBC, [Rmin,Rmax], ...
    'ObjectPolarity','bright','Sensitivity',0.85);

%% Watershed method
% tao hinh anh RBC chua nhan dien duoc
RBC_watershed1=RBC;
for k=1:size(center,1)
    R=radius(k);
    x0=center(k,1);
    y0=center(k,2);
    for i=1:size(RBC_watershed1,1)
        for j=1:size(RBC_watershed1,2)
            y=i;
            x=j;
            if (RBC_watershed1(i,j)==1)&&(sqrt((x-x0)^2+(y-y0)^2)<=R)
                RBC_watershed1(i,j)=0;
            end
        end
    end
end
%RBC_watershed1=bwareaopen(RBC_watershed1,1000); 
% Watershed
D = bwdist(~RBC_watershed1);
D = -D;
D = imgaussfilt(D, 10);
L = watershed(D);
L(~RBC_watershed1) = 0;
RBC_watershed2 = L > 0;
% tim nguong dien tich phan biet RBC voi nhieu
[B,L] = bwboundaries(RBC_watershed2,'noholes');        % tim duong bien
STATS = regionprops(RBC_watershed2, 'all');            % chi lay Area
area=zeros(size(B,1),1);
for i=1:size(B,1)
    area(i)=STATS(i).Area;
end
area=sort(area);
breakpoint=clustering(area);              % nguong cua phan nhom du lieu
areathreshold=area(breakpoint+1);         % nguong dien tich

%% Couting RBC, WBC, mask into 2-D image
s=size(center);                           % Couting RBC

n=bwconncomp(WBC);                        % Find connected components in binary image
numberWBC=n.NumObjects;                   % Couting WBC

Outline=bwperim(WBC);                     % Find perimeter of objects in binary image
Ed = imdilate(Outline,ones(4,4));         % Dilate image
WBCoverlay=imoverlay(RGB,Ed,'g');         % Burn binary mask into 2-D image

%% Imshow images and results
figure(2)
subplot(2,2,1); imshow(RGB);        title('RGB image');
subplot(2,2,2); imshow(S);          title('S-Channel');
subplot(2,2,3); imshow(RBC);        title('Segmentation RBC');
subplot(2,2,4); imshow(WBCoverlay); 

hold on
anomalous_RBC=0;
% phat hien RBC bat thuong, chong lap
for i=1:size(B,1) 
    if STATS(i).Area>=areathreshold  
        C0=(STATS(i).Perimeter)^2/STATS(i).Area;
        if (C0>=14)&&(C0<=27)
            plot(B{i}(:,2), B{i}(:,1), 'LineWidth',3,'Color', 'blue'); % duong bien
            anomalous_RBC=anomalous_RBC+1;
        end
    end
end
viscircles(center, radius,'LineStyle','-');                            % mask into 2-D image
title(strcat('RBC Hough transform:',num2str(s(1)),', RBC Watershed:',num2str(anomalous_RBC),', WBC:',num2str(numberWBC)));




