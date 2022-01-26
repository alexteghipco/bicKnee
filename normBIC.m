function [C1, Cm, C2, diffBic, iPt, iCluster, optimalCluster, optimalClusterIdx, optimalClusterBic] = normBIC(bic,nC)
% [dBic, nBic] = normBIC(bic,clusters)
%
% Accentuate the knee in a vector of BIC values by (mainly) normalizing
% the BIC value and dividing by the number of clusters. This can also be
% used with other cluster validation indices where you want to weight some
% index against the successive change in the index (i.e., weight the index
% by slope). Though note, it assumes that a lower value in the index is
% better. 
%
% Inputs: 
% bic is a N x 1 vector of BIC indices for each clustering solution
% in a series.
% nC is a N x 1 vector of the number of clusters in each of these
% clustering solution.
% if you are not using BIC, you may need to fix the way combining c1 and c2
% is handled. For instance, if your criterion values are decreasing but the
% smallest number is the best, you should set it to 'sum'. 
%
% Outputs: 
% C1 is the normalized BIC (BIC placed in range of clustering solutions)
% Cm is the ratio of normalized BIC value and the number of clusters
% C2 is Cm placed in the same range as C1 (i.e., range of clustering solutions) 
% diffBic combines C1 and C2 while accounting for multiple local maxima by
%   either adding or taking the absolute value of the differences between C1
%   and C2, depending on whether the original bic curve is increasing or
%   decreasing. These values are divided by 2 to set diffBic into the same
%   range as C1.
% iPt is the intersecting point between Normalized BIC and diffBic
% (normalized BIC will be lower than diffBic up until a certain point and
% then will mostly be greater). It excludes the lowest clustering solution.
% 
% Motivation for normalization of BIC:
% The optimal clustering solution with the BIC criterion is conventionally
% taken to be the first decisive local maximum of the BIC curve (the curve
% being comprised of BIC indices for a series of clustering solutions). In
% practice, this heuristic does not always approximate the best clustering
% solution because clustering algorithms can be inefficient or overtrained.
% The BIC index along the BIC curve that better approximates the optimal
% solution highlights a large change in BIC such that subsequent BIC
% indices reach a plateau. Thus, identifying the knee in the BIC curve is a
% better choice for deciding which local maximum has the strongest evidence
% for the correct number of clusters. However, identifying the knee can be
% challanging too because there can be several local maximums in the BIC
% curve. One approach has been to look at the slope in BIC between each
% pair of successive clustering solutions. This highlights the largest
% changes along the BIC curve but ignores the magnitude of BIC, and
% therefore the overall quality of the clustering solution. 
%
% TL;DR: We want to identify the clustering solution that drives
% the largest change along the BIC curve, but which also highlights a
% clustering solution of relatively good quality (relative to other BIC
% indices along the BIC curve). 
%
% Explanation of normalization process: To weight the BIC value by local
% changes in BIC magnitude we use the procedure outlined by Zhao et al in
% "Knee Point Detection on Bayesian Information Criterion". It involves the
% following steps: 1) Normalize the BIC values into range of clustering
% solutions which they are attached to (C1), 2) Divide by number of
% clusters (Cm), 3) normalize the resulting BIC values into the same range
% as before (i.e., range of clustering solutions; C2). Because this
% procedure seeks a large BIC value, if the original BIC curve is
% increasing this will be determined by the maximum clustering solution,
% but if it is decreasing it will be determined by the minimum clustering
% solution. Thus, step 4) is to determine which case the current BIC curve
% falls under. Step 5) is to calculate "DiffBic", which will be different
% for these two cases. In case the original BIC curve increases, the
% normalized BIC will be added to the normalized BIC weighted by the number
% of clusters (Cm) and put in the same range (C2). In the case it
% decreases, the normalized BIC will be subtracted from the weighted
% normalized BIC and the absolute value of that result will be taken.
% Finally, step 6) is to find where DiffBic intersects with normalized BIC.
%
% For the original paper, see: http://cs.uef.fi/sipu/pub/KneePointBIC-ICTAI-2008.pdf
%
% Alex Teghipco // alex.teghipco@uci.edu // 05/13/19

overrideDiff = 'auto'; % can be 'auto', 'sum', or 'diff'

%% 1) Normalize BIC values into range of clusters to get C1
for i = 1:length(bic)
    C1(i,1) = ((max(nC) - min(nC)) * (bic(i) - min(bic)))/(max(bic) - min(bic));
end

%% 2) Divide C1 by the number of clusters m to get value Cm
Cm = C1./nC; % The value of Cm calculates the ratio between the normalized BIC value and the number of clusters, which reveals the global trend of the BIC curve

%% 3) Cm is normalized into the same range of clusters as C1 was to obtain C2
for i = 1:length(Cm)
    C2(i,1) = ((max(nC) - min(nC)) * (Cm(i) - min(Cm)))/(max(Cm) - min(Cm));
end

%% ---- two cases of normalized BIC ---- %%
% The original BIC curve can have a globally increasing trend or a globally
% decreaseing trend. A large BIC value is preferred to be the optimal
% solution. In the case of a globally increasing trend, this value depends
% on the maximum cluster size, and in the case of a globally decreasing
% trend this value depends on the minimum cluster size. In the case of a
% globally increasing trend, C2 reaches several local maximums. When C2
% finds the point that indicates the largest change, it will not continue
% to increase. The largest value of C2 is considered as the most
% significant change. To reach maximum information, C1 and C2 will be
% summed. In the other case (i.e., C1 is decreasing), C2 will be decreasing
% as well. As such, the absolute subtraction of them will determine the
% most significant change. 

%% 4) Determine if Cm is globally increasing or decreasing and start plotting stuff
f = figure;
p = plot(nC,C1);

f.CurrentAxes.XTickLabelRotation = 90;
set(gcf,'color','w');
set(gca,'box','off')
p.Color = [0.7 0.7 0.7];
p.LineWidth = 2;
p.LineStyle = '--';
p.Marker = 'o';
p.MarkerSize = 7;
p.MarkerFaceColor = [1 1 1];
f.CurrentAxes.YGrid = 'on';
f.CurrentAxes.XGrid = 'on';
xlabel('Clustering solution','FontSize',18)
ylabel('BIC','FontSize',18)
f.CurrentAxes.FontSize = 12;

r = corr(bic,(1:length(bic))');

%% 5) Get DiffBic
switch overrideDiff
    case 'auto'
        if r < 0
            title(['Original BIC decreases globally according to linear correlation (r = ' num2str(r) ')'],'FontSize',28)
            diffBic = (abs(C1 - C2))./2;
        elseif r > 0
            title(['Original BIC increases globally according to linear correlation (r = ' num2str(r) ')'],'FontSize',28)
            diffBic = (C1 + C2)./2;
        end
    case 'sum'
        title(['Original BIC decreases globally according to linear correlation (r = ' num2str(r) ')'],'FontSize',28)
        diffBic = (C1 + C2)./2;
        
    case 'diff'
        title(['Original BIC increases globally according to linear correlation (r = ' num2str(r) ')'],'FontSize',28)
        diffBic = (abs(C1 - C2))./2;
end
        
%% 6) Max clustering solution refinement
% We can refine the range of clustering solutions we look at because
% everything after the knee plateaus. The original range is arbitrary
% anyway since it's difficult to produce a clustering solution for k =
% 1...maximum number of points in the data. 
hold on
p2 = plot(nC,diffBic);
p2.Color = [0.2 0.2 0.2];
p2.LineWidth = 2;
p2.LineStyle = '--';
p2.Marker = 'o';
p2.MarkerSize = 7;
p2.MarkerFaceColor = [1 1 1];
legend('Normalized BIC (C1)','Combined normalized BIC and normalized BIC weighted by number of clusters')

iPt = []; 
searchIdx = 2;
if diffBic(2) > C1(2)
    dir = 'lt'; % knee is where diffBic is less than C1
else
    dir = 'gt'; % knee is where diffBic is greater than C1
end

while searchIdx <= length(diffBic) - 1
    switch dir
        case 'lt'
            if diffBic(searchIdx) <= C1(searchIdx)
                iPt = searchIdx;
            end
        case 'gt'
            if diffBic(searchIdx) >= C1(searchIdx)
                iPt = searchIdx;
            end
    end
    if ~isempty(iPt)
        break
    end
    searchIdx = searchIdx + 1;
end

if ~isempty(iPt)
    if iPt == 1 || iPt == length(bic)
        warning('The knee point was either the first or last clustering solution. This could be wrong. Inspect figure more closely.')
    end
    [optimalClusterBic,optimalClusterIdx] = max(diffBic(1:iPt));
    optimalCluster = nC(optimalClusterIdx);
    iCluster = nC(iPt);
else
    warning('No knee detected...optimal k is set to the largest BIC value but the smallest value might make more sense')
    warning('You can try to find a subset of clustering solutions to focus on. This is especially a good idea if your BIC curve shows both increasing and decreasing global trends (i.e., where there is a clear point after which overfitting starts to occur)')
    
    [optimalClusterBic,optimalClusterIdx] = max(diffBic(1:end));
    optimalCluster = nC(optimalClusterIdx);
    iCluster = [];
end
