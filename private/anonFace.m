function anonFace(objJpg) %anonymize face by replacing skintones with grey
ogI = imread(objJpg);
I = ogI; 

% filter rgb salt and pepper
pix = 25; f = ones(pix)/pix^2;
I(:,:,1) = filter2(f,I(:,:,1));
I(:,:,2) = filter2(f,I(:,:,2));
I(:,:,3) = filter2(f,I(:,:,3));

% extract rgb as double
red  = double(I(:,:,1)); green = double(I(:,:,2)); blue = double(I(:,:,3));
% normalize by red and mask skintone values
greenR = green./red; blueR = blue./red;
gMask = (greenR>0.35 & greenR<0.85); uMask = (blueR>0.1 &  blueR<0.65);
mask = gMask & uMask;
ogI(mask(:,:,[1,1,1])) = 128;
% Overwrite Image
imwrite(ogI, objJpg)
end