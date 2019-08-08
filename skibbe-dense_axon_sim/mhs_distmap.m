function D=mhs_distmap(A,shape,varargin)


    maxdist=3;
    consider_scale=false;

    sigma=2.5;
    gamma=1;
    mode=0;
    lininterp=true;
    vscale=[1,1,1];
    weight=false;
    normalize=true;
    
    offset=[0,0,0];
    for k = 1:2:length(varargin)
            eval(sprintf('%s=varargin{k+1};',varargin{k}));
    end


    if ~exist('pathid','var')
       idx=reshape((A.connections(1:2,:)+1),[1,2*size(A.connections,2)]);
       dataP=A.data([1:3,8,21],idx);
    else
       path_ids=A.data(21,:);
       valid=(path_ids(A.connections(1,:)+1)==path_id);
       idx=reshape((A.connections(1:2,valid)+1),[1,2*sum(valid)]);
       dataP=A.data([1:3,8,21],idx);
    end;


    D=mhs_distmapC(single(dataP),{'shape',shape([3,2,1]),'normalize',normalize,'b',maxdist,'lininterp',lininterp,'gamma',gamma,'sigma',sigma,'consider_scale',consider_scale,'mode',mode,'vscale',vscale([3,2,1]),'offset',offset([3,2,1]),'weight',weight});




% 
% 
% 
% % if nargin<3
% %     maxdist=5;
% % end;
% % 
% % if nargin<4
% %     consider_scale=false;
% % end;
% 
%   maxdist=3;
%   consider_scale=false;
% 
% if ~exist('pathid','var')
%    idx=reshape((A.connections(1:2,:)+1),[1,2*size(A.connections,2)]);
%    dataP=A.data([1:3,8,21],idx);
% else
%    path_ids=A.data(21,:);
%    valid=(path_ids(A.connections(1,:)+1)==path_id);
%    idx=reshape((A.connections(1:2,valid)+1),[1,2*sum(valid)]);
%    dataP=A.data([1:3,8,21],idx);
% end;
% 
% 
% 
% if nargout>1
%     [D,C]=mhs_distmapC(single(dataP),{'shape',shape([3,2,1]),'b',maxdist,'consider_scale',consider_scale,'rgb',true});
% 
% else
%     D=mhs_distmapC(dataP,{'shape',shape([3,2,1]),'b',maxdist,'consider_scale',consider_scale,'rgb',false});
% end
