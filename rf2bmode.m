function modeB = rf2bmode(RF, increase)
%% rf2modeb Converts a RF image into a mode B one.
%  modeB = RF2modeB(RF, increase)
%
%         RF: RF image
%   increase: increase factor to adjust contrast
%
%  ----------------------
% | Code:  Renaud Morin  |
% | Email: morin@irit.fr |
% | Date:  June 2011     |
%  ----------------------
%
% Last modif : 10-01/2012

clear modeB


%% Initialization
if nargin==1
    increase=0;
elseif nargin<1 || nargin>2
    error('There must be 1 or 2 arguments.')
end


%% Computation
for i=1:size(RF,3)
    modeB_temp=20*log(abs(hilbert(RF(:,:,i)))+increase);
    modeB_temp=modeB_temp-min(modeB_temp(:));
    max_modeB=max(max(modeB_temp));
    modeB_temp = 255.*modeB_temp/max_modeB;
%     modeB_temp = uint8(255.*modeB_temp/max_modeB);
    modeB(:,:,i)=modeB_temp;
end