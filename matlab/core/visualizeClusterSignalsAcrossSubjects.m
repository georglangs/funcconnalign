function [Signal_Out] = visualizeClusterSignalsAcrossSubjects(map,BOLD,LabelName)
% function to visualize the mean BOLD signals of clusters in different
% subjects. It assumes that the cluster labels for different subjects
% correspond to each other.
%
% map{i}.(LabelName) ... cluster labels
% BOLD{i} ... Bold signals for all nodes in map i
%
% georg langs
% langs@csail.mit.edu
% 3.12.2010
%
%%
S = length(map);
for count = 1:S
C(count) = max(map{count}.(LabelName));
end
C = max(C);
%%
clf
axescounter = 1;
ha = tight_subplot(S,C,0,.03,.03);
for SubjectID = 1:S
   for ClusterID = 1:C
       if sum(map{SubjectID}.(LabelName)==ClusterID)>0
      Signal = BOLD{SubjectID}(map{SubjectID}.(LabelName)==ClusterID,:);
      Signal_Out(SubjectID,:,ClusterID) = mean(Signal);
      %subplot(S,C,axescounter)
      axes(ha(axescounter)); 
      plot(median(Signal));
       else
           Signal=1
           axes(ha(axescounter)); 
           plot(0,0)
       end
       set(gca,'XTick',[],'YTick',[]);
       set(gca,'XLim',[0 size(Signal,2)]);
       grid(gca,'on')
      axescounter = axescounter + 1;
   end
end
set(gcf,'Color','white')

