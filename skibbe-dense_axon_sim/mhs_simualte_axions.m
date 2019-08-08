% example usage:
% data=mhs_simualte_axions('shape',[512,512,512],'mode','chaos','density',1,'max_b',1,'axon_thickness',[2,3],'thick_min',1.5,'flatty_axions',true);
% data=mhs_simualte_axions('shape',[256,256,512],'mode','bundle','density',5,'dir_noise_bias',0.005,'bundle_rad',0.75);
%
% parameter, default
% 'shape',[512,512,128]  Image size
% 'min_particles', 30    discard very short candidates
% 'create_images', true  render the images
% 'axon_thickness', [2,3] thickness diameter interval for seeds
% 'thick_min',1.5  min scale
% 'density',5  as lower as more axons
% 'max_b', 1  maximum number of bifurcations per axon
% 'mode', 'chaos'  can be 'chaos' or 'bundle'
% 'shotnoise',4 photon noise, as lower, as higher
% 'background_noise',0.1 some background noise

% returns struct where
% data.GT is the graph structure with
% data.GT.data array of nodes 
%         dim 1;3 positions
%         dim 8   scale
%         dim 21 pathid
% data.GT.connections array of connections where
%         dim 1,2 node indeces (c style, for matlab +1)    
% data.Img image without noise
% data.ImgD image with noise

function [GT,Img,ImgD]=mhs_simualte_axions(varargin)

rng('shuffle') 
shape=[512,512,128];
max_d=500;
min_particles=30;

create_images=true;

axon_thickness=[2,3];
thick_min=1;   
thick_noise_fact=0.25;
var_intensity=false;

density=5;

two_dirs=true;

bif_prob=0.1;
 min_b_dist=10;
max_b=1;
 
bundle_rad=0.5;
mode='chaos';
flatty_axions=true;

background_noise=0.1;
shotnoise=4;

fprop=[0.5,1,2];

dir_noise_bias=[1,1,1];

directions=[1,1,1];

for k = 1:2:length(varargin),
            eval(sprintf('%s=varargin{k+1};',varargin{k}));
end;


if ~exist('GT','var')


 
    
switch mode
    case 'chaos'
        
        min_particles=min_particles.*shape/max(shape);
        count=1;

        
        %% seeds on all surfaces
        if directions(1)   
            for x=density:density:shape(1)-density
                
                  thick_start=max(rand*(axon_thickness(2)-axon_thickness(1)))+axon_thickness(1);

                 if (rand>0.5) || (~two_dirs)
                    ndir=[0,1,0];
                    start_pos=[x,0,rand*shape(3)];
                else
                    ndir=[0,-1,0];
                    start_pos=[x,shape(2),rand*shape(3)];
                end;
                
                if flatty_axions
                    dir_noise_bias_=dir_noise_bias;
                    dir_noise_bias_(randi(3))=0.01;
                    ndir_noise=0.01+(1-ndir).*0.1.*dir_noise_bias_;
                else
                    ndir_noise=0.01+(1-ndir).*0.1.*dir_noise_bias;
                end;
                
                
                
                A=mhs_treegen('shape',shape,'start_pos',start_pos,'ndir',ndir,'max_d',max_d,'ndir_noise',ndir_noise,'bif_prob',bif_prob,'min_b_dist',min_b_dist,'max_b',max_b,'fprop',fprop,'thick_start',thick_start,'thick_min',thick_min,'thick_noise_fact',thick_noise_fact);


                if size(A.data,2)>min_particles(1)
                    if (count==1)
                         A.data(21,:)=1;
                         C=A;
                    else
                         C=append(C,A);
                    end;
                    count=count+1;
                end;
            end;
        end;

        if directions(2)   
            for y=density:density:shape(2)-density

                thick_start=max(rand*(axon_thickness(2)-axon_thickness(1)))+axon_thickness(1);
                if (rand>0.5) || (~two_dirs)
                    ndir=[1,0,0];
                    start_pos=[0,y,rand*shape(3)];
                else
                    ndir=[-1,0,0];
                    start_pos=[shape(1),y,rand*shape(3)];
                end;

                if flatty_axions
                    dir_noise_bias_=dir_noise_bias;
                    dir_noise_bias_(randi(3))=0.01;
                    ndir_noise=0.01+(1-ndir).*0.1.*dir_noise_bias_;
                else
                    ndir_noise=0.01+(1-ndir).*0.1.*dir_noise_bias;
                end;

                A=mhs_treegen('shape',shape,'start_pos',start_pos,'ndir',ndir,'max_d',max_d,'ndir_noise',ndir_noise,'bif_prob',bif_prob,'min_b_dist',min_b_dist,'max_b',max_b,'thick_start',thick_start,'thick_min',thick_min,'thick_noise_fact',thick_noise_fact);

                if size(A.data,2)>min_particles(2)
                    if (count==1)
                         A.data(21,:)=1;
                         C=A;
                    else
                         C=append(C,A);
                    end;
                    count=count+1;
                end;
            end;
        end;


        if directions(3)   
            for x=density:density:shape(1)-density

                thick_start=max(rand*(axon_thickness(2)-axon_thickness(1)))+axon_thickness(1);
                if (rand>0.5) || (~two_dirs)
                    ndir=[0,0,1];
                    start_pos=[x,rand*shape(2),0];
                else
                    ndir=[0,0,-1];
                    start_pos=[x,rand*shape(2),shape(3)];
                end;
                if flatty_axions
                    dir_noise_bias_=dir_noise_bias;
                    dir_noise_bias_(randi(3))=0.01;
                    ndir_noise=0.01+(1-ndir).*0.1.*dir_noise_bias_;
                else
                    ndir_noise=0.01+(1-ndir).*0.1.*dir_noise_bias;
                end;
                
                
                A=mhs_treegen('shape',shape,'start_pos',start_pos,'ndir',ndir,'max_d',max_d,'ndir_noise',ndir_noise,'bif_prob',bif_prob,'min_b_dist',min_b_dist,'max_b',max_b,'thick_start',thick_start,'thick_min',thick_min,'thick_noise_fact',thick_noise_fact);

                if size(A.data,2)>min_particles(3)
                    if (count==1)
                         A.data(21,:)=1;
                         C=A;
                    else
                         C=append(C,A);
                    end;
                    count=count+1;
                end;
            end;
        end;
        
        case 'bundle'
            
        %% seeds only on one surface
        
        min_particles=min_particles.*shape/max(shape)/5;
        count=1;
        
        [X Y ] = ndgrid(0:(shape(1)-1),0:(shape(2)-1));
        mask=(mod(X,round(density))==0) & (mod(Y,round(density))==0);
        
        X2 = 2*(X - ceil(shape(1)/2))/shape(1);
        Y2 = 2*(Y - ceil(shape(2)/2))/shape(2);
        
        R2=X2.^2+Y2.^2;
        mask=(R2<bundle_rad^2) & mask;
        
        totalcand=sum((R2(:)<bundle_rad^2));
        
        
        
        
        start_poss=[X(mask),Y(mask)];
        start_poss=start_poss+randn(size(start_poss))*density*0.25;
        
        fprintf('%d startpos candidates (of %d, %d percent)\n',size(start_poss,1),totalcand,ceil(100*size(start_poss,1)/totalcand));
        
        for s=1:size(start_poss,1)
            
             thick_start=max(rand*(axon_thickness(2)-axon_thickness(1)))+axon_thickness(1);
             
            ndir=[0,0,1];
            start_pos=([start_poss(s,:),0]);
            
            
             ndir=ndir+dir_noise_bias.*randn(size(ndir));
             ndir=ndir./norm(ndir);
             
             ndir_noise=dir_noise_bias;
             max_b=0;
            
             A=mhs_treegen('shape',shape,'start_pos',start_pos,'ndir',ndir,'max_d',max_d,'ndir_noise',ndir_noise,'bif_prob',bif_prob,'min_b_dist',min_b_dist,'max_b',max_b,'thick_start',thick_start,'thick_min',thick_min,'thick_noise_fact',thick_noise_fact);

            if size(A.data,2)>min_particles(1)
                if (count==1)
                     A.data(21,:)=1;
                     C=A;
                else
                     C=append(C,A);
                end;
                count=count+1;
            end;
        end;
        
        
end;        

        GT=C;
end;        


%% tree2image
if create_images
    D=mhs_distmap(GT,shape,'maxdist',5);
    Img=single(exp(-D.^4));
    ImgD=single(mhs_poisson_noise(max(Img+background_noise+background_noise*randn(size(Img)),0.0),shotnoise));
    path_ids=unique(GT.data(21,:));
    rindx=[1:numel(path_ids)];
    
end;


if nargout==1
    data.GT=GT;
    if create_images
        data.Img=Img;
        data.ImgD=ImgD;
        if var_intensity
            data.ImgD2=ImgD2;
            data.ImgW=ImgW;
        end;
        data.rindx=rindx;
        data.params=varargin;
        GT=data;
    end;
end;

function C=append(C,B)

pid=max(C.data(21,:));


offset=size(C.data,2);

tmp=B.data;
tmp(21,:)=pid+1;
C.data=cat(2,C.data,tmp);
tmp=B.connections;
tmp(1:2,:)=tmp(1:2,:)+offset;
C.connections=cat(2,C.connections,tmp);


% function n_img=mhs_poisson_noise(A, sn)
%     scale = 1.2 * 10 ^ -3;
%     p_n = poissrnd(scale, size(A));
%     n_img = A + p_n;
function img=mhs_poisson_noise(img,lscale)
    img=poissrnd(img*lscale)/lscale;
