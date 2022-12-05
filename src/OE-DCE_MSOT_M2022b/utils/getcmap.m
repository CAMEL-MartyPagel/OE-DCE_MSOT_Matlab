function fgmap = new_getcmap(cmap)
% convert cmap string to numerical



if isempty(cmap), % display list
    fgmap = {'contrast','autumn','spring','oxy','cyan-darkblue','red','green','cyan','yellow-darkred','pink','orange','yellow','rainbow','jet','parula','hot','bone','gray','ampel','green-yellow','light yellow','light blu','light red','light green','magma','inferno','plasma','viridis','blue'}';
elseif strcmp(cmap,'gray'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,1) = lin;    
    fgmap(:,2) = lin;    
    fgmap(:,3) = lin;    
elseif strcmp(cmap,'yellow'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,1) = lin;    
    fgmap(:,2) = lin;    
% elseif strcmp(cmap,'pink'),
%    fgmap = zeros(64,3);
%    lin = linspace(0,1,64);
%    fgmap(:,1) = lin;
%    fgmap(:,3) = lin*0.5;
elseif strcmp(cmap,'autumn'),
   fgmap = zeros(64,3);
   lin = linspace(0,1,32);
    fgmap(1:32,1) = lin;
    fgmap(33:64,1) = 1;
    fgmap(33:64,2) = lin;
elseif strcmp(cmap,'contrast'),
   fgmap = zeros(64,3);
   lin = linspace(0,1,32);
    fgmap(1:32,3) = lin;
    fgmap(33:64,1) = lin;
    fgmap(33:64,2) = lin;
    fgmap(33:64,3) = flip(lin);
elseif strcmp(cmap,'green'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,2) = lin;
elseif strcmp(cmap,'red'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,1) = lin;
elseif strcmp(cmap,'blue'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,3) = lin;
elseif strcmp(cmap,'cyan-darkblue'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,1) = 0;%0.5*flip(lin);
    fgmap(:,2) = 1*flip(lin);
    fgmap(:,3) = 0.99;
%     fgmap
elseif strcmp(cmap,'cyan'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,2) = lin;
    fgmap(:,3) = lin;
elseif strcmp(cmap,'yellow-darkred'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,3) = 0*flip(lin);
    lin = linspace(0,1,32);
    fgmap(1:32,2) = flip(lin);
    fgmap(1:32,1) = 1;
    fgmap(33:64,1) = 1-0.5*lin;%0.55+0.45*flip(lin);
%     lin = (linspace(0,1,32));
%     fgmap(33:64,3) = 0.25*(lin);
%     fgmap(33:64,2) = 0.25*(lin);
%     fgmap(1:32,1) = lin;
%     fgmap(33:64,1) = 1;
%     fgmap
%     fgmap(:,1) = lin;
elseif strcmp(cmap,'orange'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,64);
    fgmap(:,2) = 0.2+0.8*flip(lin);
    fgmap(:,1) = 1;
    fgmap(:,3) = 0;

elseif strcmp(cmap,'spring'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,32);
    fgmap(33:64,2) = lin;
    fgmap(1:32,3) = lin; 
    fgmap(33:64,3) = fliplr(lin);
elseif strcmp(cmap,'oxy'),
    fgmap = ones(64,3);
    lin = (linspace(0,1,32)).^2;
    fgmap(1:32,1) = lin;
    fgmap(1:32,2) = lin; 
    fgmap(33:64,2) = fliplr(lin);
    fgmap(33:64,3) = fliplr(lin);
elseif strcmp(cmap,'rainbow'),
    fgmap = ones(64,3);
    lin = linspace(0,1,32);
    fgmap(1:32,1) = lin;
    fgmap(1:32,3) = fliplr(lin); 
    fgmap(33:64,2) = fliplr(lin);
    fgmap(33:64,3) = lin;
elseif strcmp(cmap,'ampel'),
    fgmap = zeros(64,3);
    lin = linspace(0,1,39);
    fgmap(1:39,2) = lin;
    fgmap(40:64,2) = 1;
    lin = linspace(0,1,19);
    fgmap(1:39,1) = 1;
    fgmap(40:58,1) = fliplr(lin); 
elseif strcmp(cmap,'green-yellow'),
    fgmap = buildcmap('kgy');
elseif strcmp(cmap,'light yellow'),
    fgmap = buildcmap('kyw');
    
elseif strcmp(cmap,'light yellow'),
    fgmap = buildcmap('kyw');
elseif strcmp(cmap,'light blu'),
    fgmap = buildcmap('kbw');
elseif strcmp(cmap,'light red'),
    fgmap = buildcmap('krw');
elseif strcmp(cmap,'light green'),
    fgmap = buildcmap('kgw');
    
elseif strcmp(cmap,'grey'),
    fgmap = eval('gray');    
elseif strcmp(cmap,'magma'),
    fgmap=magma(64);
elseif strcmp(cmap,'inferno'),
    fgmap=inferno(64);
elseif strcmp(cmap,'plasma'),
    fgmap=plasma(64);
elseif strcmp(cmap,'viridis'),
    fgmap=plasma(64);    
elseif ischar(cmap),
    fgmap = eval(cmap);
elseif isnumeric(cmap) && size(cmap,2) >= 3,
    fgmap = cmap;
else
    error('Unknown input for getcmap');
    fgmap = [];
end