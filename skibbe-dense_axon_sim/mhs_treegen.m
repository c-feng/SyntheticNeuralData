
function [A,t_tree]=mhs_treegen(varargin)


    shape=[32,128,128];
    max_d=25;
    min_b_dist=10;
    max_b=2;
    ndir=[0,0,1];
    ndir_noise=[0.1,0.1,0.1];
    bif_prob=0.1;
    
    bif_rule=1;
    
    thick_start=2;
    thick_min=0.5;
    thick_noise_fact=0.25;
    global_curve=0;
    step=5;
    max_depth=inf;
    
    for k = 1:2:length(varargin),
            eval(sprintf('%s=varargin{k+1};',varargin{k}));
    end;
    
    
    if ~exist('start_pos','var')
        t_root.pos=[shape(1)/2,shape(2)/2,0];    
    else
        t_root.pos=start_pos;    
    end;
    
    t_root.dir=ndir;
    t_root.scale=thick_start/2;
    
    t_root.pr={};
    t_root.su={};
    t_root.depth=0;
    
    
    t_tree{1}=t_root;
    
    data.max_d=max_d;
    data.shape=shape;    
    data.min_b_dist=min_b_dist;
    data.max_b=max_b;
    data.ndir_noise=ndir_noise;
    data.bif_prob=bif_prob;
    data.thickness=thick_start;
    data.thick_min=thick_min;
    data.thick_noise_fact=thick_noise_fact;
    data.rootdir=ndir;
    data.global_curve=global_curve;
    data.bif_rule=bif_rule;
    data.step=step;
    data.max_depth=max_depth;
    
    t_tree=gen_tree(t_tree,data);
    A=gen_ntree(t_tree);
    A.data(21,:)=1;
    
    
    function [d1,d2]=myelinated_axion_bif_dia(data)
        d0=data.thickness;
        d1=d0*(1-rand*0.75);
        v=1; 
        d2=(d0^(v+2) - d1^(v+2))^(1/(v+2));
        
    
    function [t_tree,c]=gen_tree(t_tree,data,d,b,c,bdist)
        
        
        if nargin==2
            d=0;
            b=0;
            c=1;
            bdist=0;
        end;
        
                bdist=bdist-1;
        
        if (d>=data.max_d)
            return;
        end;
        
        step=data.step;
        bif_prob=data.bif_prob*(1/(10*(b+2)));
        dir_noise=data.ndir_noise;
        
        dir_noise_bif=[0.5,0.1];
        
        if (d<2)
            bif_prob=0;
        end;
        

        current_node=t_tree{c};
        
        bif_created=false;
        
        
        
        % create bif
        if rand<bif_prob && (bdist<=0) && (b<data.max_b) 
            if (data.thickness>data.thick_min)
            
            dir_o=cross(current_node.dir,[1,0,0]);
            
            n_dir0=current_node.dir;
            n_dir0=n_dir0*(5+randn*dir_noise_bif(1))+dir_o*(1+randn*dir_noise_bif(2));
            n_dir0=n_dir0./norm(n_dir0);
            
            n_dir1=current_node.dir;
            n_dir1=n_dir1*(5+randn*dir_noise_bif(1))-dir_o*(1+randn*dir_noise_bif(2));
            n_dir1=n_dir1./norm(n_dir1);
            
            fprintf('proposing bif %f %f\n',dot(n_dir1,n_dir0),dot(n_dir1,current_node.dir));

            
            switch (data.bif_rule)

                case 1
                        [d1,d2]=myelinated_axion_bif_dia(data);

                        [pos_a,dir_a,pos_b,dir_b]=create_bif(current_node.pos,current_node.dir,data.thickness,d1,d2); 

                        data_new{1}=data;
                        data_new{2}=data;
                        data_new{1}.thickness=d1;
                        data_new{2}.thickness=d2;
                        

                        bdist=data.min_b_dist;
                        n_dirs={dir_a,dir_b};
                        n_poss={pos_a,pos_b};
                        c_back=c;
                        for e=1:2
                            n_dir=n_dirs{e};
                            n_pos=n_poss{e};
                            data_new{e}.rootdir=n_dir;
                            next_node.dir=n_dir;
                            next_node.pos=n_pos;
                            next_node.pr={c_back};
                            next_node.su={};
                            next_node.depth=current_node.depth+norm(next_node.pos-current_node.pos);

                            next_node.scale=create_scale(data_new{e});

                            if any(next_node.pos<0) || any(next_node.pos>(data.shape-1)) || ( next_node.depth>data.max_depth)
                                return;
                            end;

                            t_tree{c_back}.su=[t_tree{c_back}.su,c+1];
                            t_tree{c+1}=next_node;
                            [t_tree,c]=gen_tree(t_tree,data_new{e},d+1,b+1,c+1,bdist);
                        end;
                        fprintf('creating bif %d\n',c_back);
                        bif_created=true;
                    

                    case 2
                        [d1,d2]=myelinated_axion_bif_dia(data);


                        data_new{1}=data;
                        data_new{2}=data;
                        data_new{1}.thickness=d1;
                        data_new{2}.thickness=d2;

                        bdist=data.min_b_dist;
                        n_dirs={n_dir0,n_dir1};
                        c_back=c;
                        for e=1:2
                            n_dir=n_dirs{e};
                            data_new{e}.rootdir=n_dir;
                            next_node.dir=n_dir;
                            next_node.pos=current_node.pos+n_dir*step;
                            next_node.pr={c_back};
                            next_node.su={};
                            next_node.depth=current_node.depth+norm(next_node.pos-current_node.pos);

                            next_node.scale=create_scale(data_new{e});

                            if any(next_node.pos<0) || any(next_node.pos>(data.shape-1)) || ( next_node.depth>data.max_depth)
                                return;
                            end;

                            t_tree{c_back}.su=[t_tree{c_back}.su,c+1];
                            t_tree{c+1}=next_node;
                            [t_tree,c]=gen_tree(t_tree,data_new{e},d+1,b+1,c+1,bdist);
                        end;
                        fprintf('creating bif %d\n',c_back);
                        bif_created=true;
                    
            end;
            else
                fprintf('too thin, cannot create bif\n');
            end;
        end;
        
        % just adding a segment    
        if ~bif_created
            n_dir=current_node.dir;
            
            
            
            n_dir_proposed=[];
            valid_dir=false;
            while isempty(n_dir_proposed) || ~(valid_dir)
                n_dir_proposed=n_dir+dir_noise.*randn(size(n_dir));
                n_dir_proposed=n_dir_proposed./norm(n_dir_proposed);
                
                if dot(n_dir_proposed,data.rootdir)>data.global_curve
                    valid_dir=true;
                end;
            end;
            n_dir=n_dir_proposed;
            
            next_node.dir=n_dir;
            next_node.pos=current_node.pos+n_dir*step;
            next_node.pr={c};
            next_node.su={};
            next_node.depth=current_node.depth+1;
            next_node.scale=create_scale(data);
            next_node.depth=current_node.depth+norm(next_node.pos-current_node.pos);

            if any(next_node.pos<0) || any(next_node.pos>(data.shape-1)) || ( next_node.depth>data.max_depth)
                return;
            end;
            
            
            
            t_tree{c}.su=[current_node.su,c+1];
            t_tree{c+1}=next_node;
            [t_tree,c]=gen_tree(t_tree,data,d+1,b,c+1,bdist);
        end;
        
        

        

function s=create_scale(data)        
    
    s=data.thickness+data.thick_noise_fact*data.thickness*randn;
    s=max(s,data.thick_min)/2;
    
    return;
    
    s=s+data.s*randn;
    s=max(s,opts.min);
    s=min(s,opts.max);
    
    return
    s=s+0.5*randn;
    s=max(s,1);
    s=min(s,2);
    
        

function  [pos_a,dir_a,pos_b,dir_b]=create_bif(pos_c,dir_c,c,a,b)   
        [vn1,vn2]=createOrths(dir_c);
        
        P1=[c,0];
        w=(b^2+c^2-a^2)/(2*b*c);
        P2=[b*w,b*sqrt(1-w^2)];
    
        
        
        Pc=P1/2;
        
        Pa=(P1+P2)/2-Pc;
        Pb=P2/2-Pc;
        Pc=Pc-Pc;
        
        Pn=cross([P2,0],[P1,0]);
        Pn=Pn/norm(Pn);
        
        Pnb=cross([P2,0],Pn);
        Pna=cross([P2-P1,0],Pn);
        Pnb=Pnb/norm(Pnb);
        Pna=Pna/norm(Pna);
        
        dir_co=randn*vn1+randn*vn2;
        dir_co=dir_co/norm(dir_co);
        
        offset=1+a/2;
        pos_a=offset*(Pa(1)*dir_co+Pa(2)*dir_c)+pos_c;
        offset=1+b/2;
        pos_b=offset*(Pb(1)*dir_co+Pb(2)*dir_c)+pos_c;
        
        dir_a=Pna(1)*dir_co+Pna(2)*dir_c;
        dir_b=Pnb(1)*dir_co+Pnb(2)*dir_c;
    
        if dot(dir_a,dir_c)<0
            dir_a=-dir_a;
        end;
        if dot(dir_b,dir_c)<0
            dir_b=-dir_b;
        end;
    
function [vn1,vn2]=createOrths(v)
        if ((v(1))>abs(v(2)))
            if (abs(v(2))>abs(v(3)))
              vn1=cross(v,[0,0,1]);
            else
              vn1=cross(v,[0,1,0]);
            end;
        else
            if (abs(v(1))>abs(v(3)))
              vn1=cross(v,[0,0,1]);
            else
              vn1=cross(v,[1,0,0]);
            end;      
        end;
      
      vn1=vn1/norm(vn1);
      vn2=cross(vn1,v);
    
    
        
function A=gen_ntree(t_tree)
             
        npts=numel(t_tree);    
            
        A.data=zeros(25,npts);
        
        for p=1:npts
            A.data(1:3,p)=t_tree{p}.pos;
            A.data(4:6,p)=t_tree{p}.dir;
            A.data(8,p)=t_tree{p}.scale;
        end;
        
        nedges=0;
        for p=1:npts
            nedges=nedges+numel(t_tree{p}.su);
        end;
        
        A.connections=zeros(6,nedges);
        count=1;
        for p=1:npts
            for c=1:numel(t_tree{p}.su);
                ptA=p;
                ptB=t_tree{p}.su{c};
                A.connections(1:2,count)=[ptA,ptB]-1;
                A.data(22,A.connections(2,count)+1)=1;
                A.data(22,A.connections(1,count)+1)=1;
                count=count+1;
            end;
        end;
        
    
    