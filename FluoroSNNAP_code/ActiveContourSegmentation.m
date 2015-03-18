function L = ActiveContourSegmentation(I)
% Segmentation by active-contour evolution
% User selects a seed point with a mouse click
% Use the seed point to generate a mask and perform activecontour evolution
% to segment the object (neuron cell body). Display the cell boundary
%%
I=uint16(I);
figure; imshow(imadjust(I));
L = false(size(I));
hold on
while(1)
    h = impoint;
    mask=false(size(I));
    coord = floor(h.getPosition);
    mask(coord(2)-4:coord(2)+4,coord(1)-4:coord(1)+4)=1;
    bw = activecontour(I,mask);
    L = L + bw;
    delete(h);
    B = bwboundaries(bw);
    for i=1:length(B)
        plot(B{i}(:,2),B{i}(:,1),'r');
    end
end