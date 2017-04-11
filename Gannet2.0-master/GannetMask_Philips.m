function [MRS_struct ] = GannetMask_Philips(sparname, nii_file, MRS_struct, ii)


if (nargin == 2)
    MRS_struct.ii=1;
    ii = 1;
end


% this relies on SPM, nifti exported by Philips, and spar/sdat

% some code to make ok wth 2 or 3  inputs
% fix nifti inputs and writing output mask to directory of where the nifti
% is. 

[pathspar,namespar,ext] = fileparts(sparname);
[pathnii,namenii,extnii] = fileparts(nii_file);

fidoutmask = fullfile(pathnii,[namespar '_mask.nii'])



sparheadinfo = textread(sparname, '%s');

sparidx=find(ismember(sparheadinfo, 'ap_size')==1);
MRS_struct.p.voxsize(ii,2) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'lr_size')==1);
MRS_struct.p.voxsize(ii,1) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'cc_size')==1);
MRS_struct.p.voxsize(ii,3) = str2num(sparheadinfo{sparidx+2});

sparidx=find(ismember(sparheadinfo, 'ap_off_center')==1);
MRS_struct.p.voxoff(ii,2) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'lr_off_center')==1);
MRS_struct.p.voxoff(ii,1) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'cc_off_center')==1);
MRS_struct.p.voxoff(ii,3) = str2num(sparheadinfo{sparidx+2});

sparidx=find(ismember(sparheadinfo, 'ap_angulation')==1);
MRS_struct.p.voxang(ii,2) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'lr_angulation')==1);
MRS_struct.p.voxang(ii,1) = str2num(sparheadinfo{sparidx+2});
sparidx=find(ismember(sparheadinfo, 'cc_angulation')==1);
MRS_struct.p.voxang(ii,3) = str2num(sparheadinfo{sparidx+2});


V=spm_vol(nii_file);
[T1,XYZ]=spm_read_vols(V);
H=spm_read_hdr(nii_file);
%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
halfpixshift = -H.dime.pixdim(1:3).'/2;
halfpixshift(3) = -halfpixshift(3);
XYZ=XYZ+repmat(halfpixshift,[1 size(XYZ,2)]);

% get information from SPAR - change later to be read in
% 
ap_size = MRS_struct.p.voxsize(ii, 2);
lr_size = MRS_struct.p.voxsize(ii, 1);
cc_size = MRS_struct.p.voxsize(ii, 3);
ap_off = MRS_struct.p.voxoff(ii, 2);
lr_off = MRS_struct.p.voxoff(ii, 1);
cc_off = MRS_struct.p.voxoff(ii, 3);
ap_ang = MRS_struct.p.voxang(ii, 2);
lr_ang = MRS_struct.p.voxang(ii, 1);
cc_ang = MRS_struct.p.voxang(ii, 3);
% 
% 
%We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;


ap_ang = -ap_ang;
lr_ang = -lr_ang;



% define the voxel - use x y z  
% currently have spar convention that have in AUD voxel - will need to
% check for everything in future...
% x - left = positive
% y - posterior = postive
% z - superior = positive
vox_ctr = ...
      [lr_size/2 -ap_size/2 cc_size/2 ;
       -lr_size/2 -ap_size/2 cc_size/2 ;
       -lr_size/2 ap_size/2 cc_size/2 ;
       lr_size/2 ap_size/2 cc_size/2 ;
       -lr_size/2 ap_size/2 -cc_size/2 ;
       lr_size/2 ap_size/2 -cc_size/2 ;
       lr_size/2 -ap_size/2 -cc_size/2 ;
       -lr_size/2 -ap_size/2 -cc_size/2 ];
   
% make rotations on voxel
rad = pi/180;
initrot = zeros(3,3);

xrot = initrot;
xrot(1,1) = 1;
xrot(2,2) = cos(lr_ang *rad);
xrot(2,3) =-sin(lr_ang*rad);
xrot(3,2) = sin(lr_ang*rad);
xrot(3,3) = cos(lr_ang*rad);

yrot = initrot;
yrot(1,1) = cos(ap_ang*rad);
yrot(1,3) = sin(ap_ang*rad);
yrot(2,2) = 1;
yrot(3,1) = -sin(ap_ang*rad);
yrot(3,3) = cos(ap_ang*rad);

zrot = initrot;
zrot(1,1) = cos(cc_ang*rad);
zrot(1,2) = -sin(cc_ang*rad);
zrot(2,1) = sin(cc_ang*rad);
zrot(2,2) = cos(cc_ang*rad);
zrot(3,3) = 1;

% rotate voxel
vox_rot = xrot*yrot*zrot*vox_ctr.';

% calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot+vox_ctr_coor;



%%% make new comment

mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
distance2voxctr=sqrt(sum((XYZ-repmat([lr_off ap_off cc_off].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr<=sphere_radius)=1;
%sphere_mask2=ones(1,(sum(sphere_mask)));

mask(sphere_mask==1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);

tri = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

mask = reshape(mask, V.dim);

V_mask.fname=[fidoutmask];
V_mask.descrip='MRS_Voxel_Mask';
V_mask.dim=V.dim;
V_mask.dt=V.dt;
V_mask.mat=V.mat;

V_mask=spm_write_vol(V_mask,mask);

T1img = T1/max(T1(:));
T1img_mas = T1img + .2*mask;

% construct output
% 
 voxel_ctr = [-lr_off -ap_off cc_off];


fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.outfile(ii,:)=fidoutmask;
%This assumes 1-mm iso T1 - need to fix at a later date.

voxel_ctr(1:2)=-voxel_ctr(1:2);

voxel_search=(XYZ(:,:)-repmat(voxel_ctr.',[1 size(XYZ,2)])).^2;
voxel_search=sqrt(sum(voxel_search,1));
[min2,index1]=min(voxel_search);

[slice(1) slice(2) slice(3)]=ind2sub( V.dim,index1);

%slice = [round(V.dim(1)/2+voxel_ctr(1)) 
 %       round(V.dim(2)/2+voxel_ctr(2)) 
 %       round(V.dim(3)/2+voxel_ctr(3))];

    
    
    
size_max=max(size(T1img_mas));
three_plane_img=zeros([size_max 3*size_max]);
im1 = squeeze(T1img_mas(:,:,slice(3)));
im1 = im1(end:-1:1,:)';  %not sure if need this '
im3 = squeeze(T1img_mas(:,slice(2),:));
im3 = im3(end:-1:1,end:-1:1)'; %may not need '
im2 = squeeze(T1img_mas(slice(1),:,:));
im2 = im2(:,end:-1:1)';

three_plane_img(:,1:size_max) = image_center(im1, size_max);
three_plane_img(:,size_max*2+(1:size_max))=image_center(im3,size_max);
three_plane_img(:,size_max+(1:size_max))=image_center(im2,size_max);

MRS_struct.mask.img(ii,:,:)=three_plane_img;
MRS_struct.mask.T1image(ii,:) = {nii_file};

%if(nargin==2)
figure(198);
imagesc(three_plane_img);
colormap('gray');
caxis([0 1])
axis equal;
axis tight;
axis off;
%end

end

