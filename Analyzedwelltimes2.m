function [output_args] = Analyzedwelltimes2(fname,framerate,numFrames)%,shutteropen)
%enter frame rate in 'sec per frame'
%AnalyzeDwlltimes2 finds time protein spends on and off Nuc. Outputs these as
%onlist and offlist, respectively. 
%The purpose of this file is to make a list of the amount of time
%spent in high and low state from idealized traces outputted from vbFRET 



file=load(fname)
disp("analyzing")
%disp(a);
%assignin('base','a',a);
FRETpath = file.path;
numPaths = size(FRETpath,2) %how many fluorophores?
%numFrames = 1200; %how many frames in each movie
statematrix = zeros(numFrames,numPaths);
FRETmatrix = zeros(numFrames,numPaths);
FRETdiff = zeros(numPaths,1);
%statestructure = struct(5,6);
for k = 1:numPaths; 
%disp(k)
a=FRETpath(1,k); 
%assignin('base','a',a);
b=unique(a{1,1});
tracelength = size(a{1,1});
%something like for legnth of trace count number in high state and low
%state
    highstate = max(b);
    lowstate = min(b);
    FRETdiff(k) = highstate-lowstate;
    assignin('base','FRETdiff',FRETdiff);
    
    
    for i = 1:tracelength
        
        if a{1,1}(i) == highstate && a{1,1}(i) ~= lowstate;
            statematrix(i,k) = 2;
            FRETmatrix(i,k) = highstate;
           % disp('highstate!')
        end
         if a{1,1}(i) == lowstate && a{1,1}(i) ~= highstate;
        
            statematrix(i,k) = 1;
            FRETmatrix(i,k) = lowstate;
            %disp('lowstate')
           
        assignin('base','statematrix',statematrix);
        assignin('base','FRETmatrix',FRETmatrix);
        %statematrix
         end
         
         if a{1,1}(i) == lowstate && a{1,1}(i) == highstate;
        
            statematrix(i,k) = 3;
            FRETmatrix(i,k) = lowstate;
            %disp('lowstate')
           
        assignin('base','statematrix',statematrix);
        assignin('base','FRETmatrix',FRETmatrix);
        %statematrix
         end
         %if there's only one state
         
assignin('base','statematrix',statematrix);
    end
end

%Now take the statelist matrix and turn it into a list. 
totcells = numFrames*numPaths;
statelist = reshape(statematrix,totcells,1);
statelist = statelist(find(statelist~=0)); 
assignin('base','statelist',statelist);
FRETlist = reshape(FRETmatrix,totcells,1);
FRETlist = FRETlist(find(FRETlist~=0)); 
assignin('base','FRETlist',FRETlist);

s.state = statelist;
s.FRET = FRETlist;
assignin('base','s',s)

numFrames = size(s.FRET,1);
lowlist = zeros(numFrames,1);
highlist = zeros(numFrames,1);
midlist = zeros(numFrames,1);

for i = 1:numFrames;
if s.state(i)  == 1;
    lowlist(i) = s.FRET(i);
end
if s.state(i) ==2;
    highlist(i) = s.FRET(i);
end

if s.state(i) == 3;
    midlist(i) = s.FRET(i);
    
end

end
%remove zeros from highlist and lowlist
highlist = highlist(find(highlist~=0)); 
lowlist = lowlist(find(lowlist~=0));
midlist = midlist(find(midlist~=0));
assignin('base','highlist',highlist);
assignin('base','lowlist',lowlist);
assignin('base','midlist',midlist);
histogram(lowlist,20);
hold on;
histogram(highlist,20);
hold on;
histogram(midlist,20);


%now analyze dwelltimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% find all transition events %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define bind and release variables
bind = zeros(3*numPaths,1);
release = zeros(3*numPaths,1);
numbind = zeros(3*numPaths,1);
numrelease = zeros(3*numPaths,1);
stitch_state = zeros(3*numPaths,2);
%stitch_state is a variable I'm trying out to try to stitch together long
%traces to get more information from them
m=1;
n=1;
for k = 1:numPaths;
    %disp(k);
m=1;
n=1;
disp('path length');
tracesize = size(FRETpath{k},1)
a=FRETpath(1,k); 
%assignin('base','a',a);
b=unique(a{1,1});
%tracelength = size(a{1,1});
%something like for legnth of trace count number in high state and low
%state
    highstate = max(b)
    lowstate = min(b)

    for i=1:(tracesize-1);%find binding events
    %disp('analyzing')
    %disp(i);
    %disp(tracelength);
    
    
  if ismember(a{1,1}(i),highstate) && ismember(a{1,1}(i+1),lowstate);
       %disp('into for loop')
       bind(k,m) = i;
       m = m+1;
       %disp('bind!!!!')
       %disp(i)
       
       %save(filename, 'bind', '-ascii', '-double', '-append')
       assignin('base','bind', bind);
    end
    %}
end
assignin('base','bind',bind);



for ii = 1:(tracesize-1)
    if ismember(a{1,1}(ii),lowstate)&&ismember(a{1,1}(ii+1),highstate);
        release(k,n) = ii;
        n=n+1;
        %disp('release!!!!')
    %save(filename, 'bind', '-ascii', '-double', '-append')
    assignin('base','bind',bind);
    assignin('base','release', release);
    end  
end
assignin('base','bind',bind)
%save('bind_values.mat','bind','-append')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% find number of bind/release events for each trace %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bindrow1 = bind(k,:);
release1 = release(k,:);
bind_row = bindrow1(find(bindrow1~=0));
bind_release = release1(find(release1~=0));
numbind(k) = length(bind_row);
numrelease(k) = length(bind_release);
assignin('base','numbind',numbind);
%assignin('base', 'lenrelease', lenrelease); 

lenbind = length(bind_row);
lenrelease = length(bind_release);

%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%find dwell times%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%case1: bind<release-->FRET: start high, bind(max)>release --> end low %%%%%%%%%%%%%%%%
%%%checked 11/6/17 I think this is correct%%%%
t=1;
r=1;
%case1:
if bind(k,1)<release(k,1)&&max(bind(k,:))>max(release(k,:))&&bind(k,1)~=0&&release(k,1)~=0;
    stitch_state(k,1) = 2; %trace starts unbound
    stitch_state(k,2) = 1; %trace ends bound
    assignin('base','stitch_state',stitch_state);
    %disp('starts unbound, ends bound');
    x=['num of binding events = ', num2str(length(bind_row))];
    y=['num of release events = ',num2str(length(bind_release))];
    disp(x);
    disp(y);
    z = ['(lenbind,lenrelease)', num2str(lenbind),num2str(lenrelease)];
    disp(z);
    for i = 1:(lenrelease);
        onmatrix(k,t) = release(k,i)-bind(k,i);
        t=t+1;
        assignin('base','onmatrix',onmatrix);
    end
    for i = 2:(lenbind);
        offmatrix(k,t) = bind(k,i) - release(k,i-1);
        r=r+1;
        assignin('base','offmatrix',offmatrix);
    end
end
%bindvalues = bind(k,:);
%releasevalues = release(k,:);

%%% starts low (first transition is up), Ends low b/c last transition is
%%% binding checked 11/6 think this is also correct
r=1;
t=1;
%case2:
if bind(k,1)>release(k,1)&&max(bind(k,:))>max(release(k,:))&&bind(k,1)~=0&&release(k,1)~=0;%&&release(k,1)~=0;
    stitch_state(k,1) = 1; %trace starts bound
    stitch_state(k,2) = 1; %trace ends bound
    %disp('starts bound, ends bound');
    x=['num of binding events = ', num2str(length(bind_row))];
    y=['num of release events = ',num2str(length(bind_release))];
    disp(x);
    disp(y);
    for i = 1:(lenbind-1);
        onmatrix(k,t) = release(k,(i+1))-bind(k,i);
        t=t+1;
        assignin('base','onmatrix',onmatrix);
    end
    for i = 1:(lenbind);
        offmatrix(k,r) = bind(k,i)-release(k,i);
        r=r+1;
        assignin('base','offmatrix',offmatrix);
    end
end
%%%%start high (first event is binding event), end high (last event is a
%%%%release)  %%10/24 good
t=1;
r=1;
%case3:
if bind(k,1)<release(k,1)&&max(bind(k,:))<max(release(k,:))&&bind(k,1)~=0&&release(k,1)~=0;
    stitch_state(k,1) = 2; %trace starts unbound
    stitch_state(k,2) = 2; %trace ends unbound
   % disp('starts unbound, end unbound');
    x=['num of binding events = ', num2str(length(bind_row))];
    y=['num of release events = ',num2str(length(bind_release))];
    disp(x);
    disp(y);
    
    for i = 1:(lenbind);
        onmatrix(k,t) = release(k,i)-bind(k,i);
        t=t+1;
        assignin('base','onmatrix',onmatrix);
    end
    
    for i = 2:(lenbind)
        %disp(j); disp(i);
        offmatrix(k,r) = bind(k,i)-release(k,i-1);
        r=r+1;
        assignin('base','offmatrix',offmatrix);
    end
    

end
%start high (release before bind), ends low (last release is after last bind)%% 
t=1;
r=1;
%case4:
if bind(k,1)>release(k,1)&&max(bind(k,:))<max(release(k,:))&&bind(k,1)~=0&&release(k,1)~=0;
    stitch_state(k,1) = 1; %trace starts bound
    stitch_state(k,2) = 2; %trace ends unbound
   % disp('starts unbound, ends bound');
    x=['num of binding events = ', num2str(length(bind_row))];
    y=['num of release events = ',num2str(length(bind_release))];
    disp(x);
    disp(y);
    
    for i = 1:(lenbind);
        onmatrix(k,t) = release(k,(i+1))-bind(k,i);
        t=t+1;
        assignin('base','onmatrix',onmatrix);
    end
    for i = 1:(lenbind);
        offmatrix(k,r) = bind(k,i)-release(k,i);
        r=r+1;
        assignin('base','offmatrix',offmatrix);
    end

end 
    
%Case5: if in low state the entire time.

%Case 6: if in high state the entire time. 
%assignin('base','bind',bind)
%}
end%end k for loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Look at spots without DNA, do I see fluctuations? %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%the hard thing about this is that I only have intensities from 
%maketracesdual2 which uses a textfile of points generated from imageJ. 
%%all other points are ignored, unfortunately. I would need to run a bunch
%%of random points through MakeTracesDual2 first. I guess try this?
%%% with this method, need to make a text file of points 'fnamerandom.txt'
%{
fname2 = strrep(fname, '(sel).mat', 'random.tif'); 
MakeTracesDual2(fname2);
%}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% find size of matricies %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[lOfM,wOfM] =size(offmatrix);
[lOnM,wOnM] =size(onmatrix);

totOnM = lOnM*wOnM;
totOfM = lOfM*wOfM;

onlist = reshape(onmatrix,totOnM,1);
offlist = reshape(offmatrix,totOfM,1);

onlist = onlist(find(onlist~=0)); 
offlist = offlist(find(offlist~=0));
assignin('base','onlist',onlist);
assignin('base','offlist',offlist);

%%now multiply the onlist and offlist by the length of the video to get the
%%dwelltime and unbound time in seconds. Framerate is frames/second
onlist = onlist*framerate;
offlist = offlist*framerate;
assignin('base','onlist',onlist);
assignin('base','offlist',offlist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% count number of blinks ('fluctuations') %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i = 1:numfluor
blinklist(i,1) = sum(offmatrix(i,:)~=0)
blinklist(i,2) = sum(onmatrix(i,:)~=0)
%blinklist(i,3) = sum(TFrelease(i,:)~=0)
i = i+1
assignin('base','blinklist', blinklist);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% plot dwell times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
histfit(onlist,30,'exponential');
title('Dwell Times');
%TextBox = uicontrol('style','text');
%set(TextBox, 'String', 'Hello World');
subplot(2,1,2)
histfit(offlist,30,'exponential');
title('Off times')

koffdist = fitdist(onlist,'exponential')
kondist = fitdist(offlist,'exponential')

%%%
%turn rate histograms into CSDs. 
koff = 1/koffdist.mu;
kon = 1/kondist.mu;

save rates.mat koff kon

xoff = linspace(0,max(onlist));
yoff = 1-exp(-koff*xoff);
koff_fit = zeros(length(xoff),2);
koff_fit(:,1) = xoff;
koff_fit(:,2) = yoff;
assignin('base','koff_fit',koff_fit);

xon = linspace(0,max(offlist));
yon = 1-exp(-kon*xon);
kon_fit = zeros(length(xon),2);
kon_fit(:,1) = xon;
kon_fit(:,2) = yon;
assignin('base','koff_fit',kon_fit);




%{
figure 
plot(xoff,yoff)
%}
dwell_list = sort(onlist);
unbound_list = sort(offlist);
dwell_list_CSD = zeros(max(onlist),2);
unbound_list_CSD = zeros(max(offlist),2);
assignin('base', 'dwell_list',dwell_list);
assignin('base', 'unbound_list',unbound_list);
tab_dwell = tabulate(dwell_list);
tab_unbound = tabulate(unbound_list);
assignin('base','tab_dwell',tab_dwell);
assignin('base','tab_unbound',tab_unbound);
disp('csd part');

%Cumsum dwell times
cum_dwell = zeros(max(onlist),2);
cum_dwell(:,1) = tab_dwell(:,1);
cum_dwell(:,2) = cumsum(tab_dwell(:,2));
%normalize
cum_dwell(:,2) = cum_dwell(:,2)/(length(dwell_list))
assignin('base','cum_dwell',cum_dwell);
%Cumsum unbound times
cum_unbound = zeros(max(offlist),2);
cum_unbound(:,1) = tab_unbound(:,1);
cum_unbound(:,2) = cumsum(tab_unbound(:,2));
%normalize
cum_unbound(:,2) = cum_unbound(:,2)/(length(unbound_list));
assignin('base','cum_unbound',cum_unbound);


figure
subplot(2,1,1) 
plot(xoff,yoff,'LineWidth', 4);
title('k_{off}');
hold on;
plot(cum_dwell(:,1),cum_dwell(:,2),'--rx','MarkerSize',.1);
hold off; 

subplot(2,1,2) 
plot(xon,yon,'LineWidth', 4);
title('k_{on}');
hold on;
plot(cum_unbound(:,1),cum_unbound(:,2),'--rx','MarkerSize',.1);
hold off; 

%now I need to find someway to output this data so that I can plot it all
%together

save koff_fit.mat koff_fit cum_dwell;

save kon_fit.mat kon_fit cum_unbound;

save 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now what I want to do is to figure out a way to stitch together traces.
%This,potentially, will give me the ability to look at dwelltimes longer
%than the life of the fluorophore or the movie. Because right now I have
%tons of binding events that extend past the length of the trace. My plan
%is to stitch together the path variables that is outputted from vbFRET. 





end %end function

