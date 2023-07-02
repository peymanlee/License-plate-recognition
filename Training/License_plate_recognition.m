% License Plate Recognition

close all
clear all
clc
load imgfildata;

%% Getting input and start processing image

img0 = imread('Images/P9170012.JPG');
figure,imshow(img0)
img = img0;
Input = rgb2gray(im2double(img));
Input = histeq(Input);
% figure,imshow(Input)
Input = adapthisteq(Input,'clipLimit',0.02,'Distribution','rayleigh');
% figure,imshow(Input)
Input = imsharpen(Input);
% Input = imadjust(Input);

threshold = graythresh(Input);
BW = imbinarize(Input,threshold);
figure
imshowpair(Input,BW,"montage")


Icorrected = imtophat(Input,strel("disk",15));
Icorrected = double(BW);
BW1 = imbinarize(Icorrected);
marker = imerode(Icorrected,strel("line",10,0));
Iclean = imreconstruct(marker,Icorrected);
Ibinary = imbinarize(Iclean);
BW2 = imcomplement(Ibinary);

%% find all rectangle shape objects in the picture

cc = bwconncomp(Ibinary);
stats = regionprops(cc, ["BoundingBox","Area"]);
roi = vertcat(stats(:).BoundingBox);
area = vertcat(stats(:).Area);
img = insertObjectAnnotation(Input,"rectangle",roi,area,"LineWidth",3);

areaConstraint = area > 1000 & area < 50000;
roi = double(roi(areaConstraint,:));
img = insertShape(Input,"rectangle",roi);

% select special size rectangles
I = Input;
width  = roi(:,3);
height = roi(:,4);
aspectRatio0 = height ./ width;
filterIdx0 = aspectRatio0 > 0.13 & aspectRatio0 < 0.5;
roi = roi( filterIdx0 ,:);
img = insertShape(I,"rectangle",roi);
% figure;imshow(img);


% remove very small or very large detected objects
picture = bwareaopen(~BW,30);
[~,cc]=size(img);
if cc>2000
    picture1=bwareaopen(picture,3500);
else
    picture1=bwareaopen(picture,3000);
end
% figure,imshow(picture1)
picture2=picture-picture1;
% figure,imshow(picture2)
picture2=bwareaopen(picture2,50);
% figure,imshow(picture2)


%% look for all words in all rectangles

% first approach
[L,Ne]=bwlabel(picture2);
% figure,imshow(picture2)
propied=regionprops(L,'BoundingBox');
bboxes0 = vertcat(propied.BoundingBox);
w0 = bboxes0(:,3);
h0 = bboxes0(:,4);
aspectRatio0 = w0./h0;
filterIdx0 = aspectRatio0' > 3;

filt0 = bboxes0(:,1) > size(I,2)*0.1 & bboxes0(:,1)+bboxes0(:,3) < size(I,2)*0.9 & bboxes0(:,2) > size(I,1)*0.1 & bboxes0(:,2)+bboxes0(:,4) < size(I,1)*0.9;
bboxes0 = bboxes0(filt0,:);


% second approach
I = im2gray(img);
I = imsharpen(I);
I = picture2;
% Detect MSER regions.
[mserRegions, mserConnComp] = detectMSERFeatures(I, ...
    "RegionAreaRange",[200 8000],"ThresholdDelta",4);
mserStats = regionprops(mserConnComp, "BoundingBox", "Eccentricity", ...
    "Solidity", "Extent", "Euler", "Image");

bbox = vertcat(mserStats.BoundingBox);
w = bbox(:,3);
h = bbox(:,4);
aspectRatio = w./h;

filterIdx = aspectRatio' > 3;
filterIdx = filterIdx | [mserStats.Eccentricity] > .995 ;
filterIdx = filterIdx | [mserStats.Solidity] < .3;
filterIdx = filterIdx | [mserStats.Extent] < 0.2 | [mserStats.Extent] > 0.9;
filterIdx = filterIdx | [mserStats.EulerNumber] < -4;

% Remove regions
mserStats(filterIdx) = [];
mserRegions(filterIdx) = [];

bboxes = vertcat(mserStats.BoundingBox);

if(length(bboxes) ~= 0)
    filt = bboxes(:,1) > size(I,2)*0.1 & bboxes(:,1)+bboxes(:,3) < size(I,2)*0.9 & bboxes(:,2) > size(I,1)*0.1 & bboxes(:,2)+bboxes(:,4) < size(I,1)*0.9;
    bboxes = bboxes(filt,:);
end


%% detecting box of licence and the alphabets

final_bbox = vertcat(bboxes0,bboxes);

max2 = 0;
sum2 = 0;
ind_i_max = 0;
for i = 1:length(roi(:,1))
    for j = 1:length(final_bbox(:,1))
        if ((final_bbox(j,1)+final_bbox(j,3)/2)>roi(i,1) && (final_bbox(j,1)+final_bbox(j,3)/2)<roi(i,1)+roi(i,3) && final_bbox(j,1)+3>roi(i,1) && final_bbox(j,1)+final_bbox(j,3)-3<roi(i,1)+roi(i,3))
            if((final_bbox(j,2)+final_bbox(j,4)/2)>roi(i,2)+roi(i,4)*0.4 && (final_bbox(j,2)+final_bbox(j,4)/2)<roi(i,2)+roi(i,4)*0.6&& final_bbox(j,2)+3>roi(i,2) && final_bbox(j,2)+final_bbox(j,4)-3<roi(i,2)+roi(i,4))
                sum2 = sum2 + 1;
            end
        end
    end
    if(sum2 > max2)
        max2 = sum2;
        ind_i_max = i;
    end
    sum2 = 0;
end

if ind_i_max ~= 0
    imgg = insertShape(img0,"rectangle",roi(ind_i_max,:));
    figure;
    imshow(imgg);

    I2 = imcrop(img0,roi(ind_i_max,:));
    figure
    imshow(I2)
    I3 = imcrop(picture2,roi(ind_i_max,:));
    figure
    imshow(I3)
    [L2,Ne2]=bwlabel(I3);
    propied2=regionprops(L2,'BoundingBox');
    hold on
    for n=1:size(propied2,1)
        rectangle('Position',propied2(n).BoundingBox,'EdgeColor','g','LineWidth',2)
    end
    hold off

    results = ocr(I3);
    txt = results.Text;
    licence_plate = regexprep(txt, '[^A-Z0-9]', '');
    

    if(length(licence_plate) == 0)
    final_output=[];
    t=[];
    for n=1:Ne2
        [r,c] = find(L2==n);
        n1=I3(min(r):max(r),min(c):max(c));
        n1=imresize(n1,[42,24]);
        x=[ ];
        totalLetters=size(imgfile,2);
        for k=1:totalLetters
            y=corr2(imgfile{1,k},n1);
            x=[x y];
        end
        t=[t max(x)];
        if max(x)>.45
            z=find(x==max(x));
            out=cell2mat(imgfile(2,z));

            final_output=[final_output out];
        end
    end

    licence_plate = final_output;
    end
    licence_plate

else
    disp('Could not find the licence plate!');
end
