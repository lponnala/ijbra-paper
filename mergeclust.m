
function [sLambda,sStart,sStop,sCz,sNz]=mergeclust(Lambda,Start,Stop,Cz,Nz)

sLambda=[]; sStart=[]; sStop=[]; sCz=[]; sNz=[];

k=1;
sLambda=[sLambda; Lambda(k)]; sStart=[sStart; Start(k)]; sStop=[sStop; Stop(k)]; sCz=[sCz; Cz(k)]; sNz=[sNz; Nz(k)];

% Lean towards higher lambda, (and if same lambda, then) smaller size, (and
% if same size, then) higher density
k=k+1;
while k<=length(Lambda)
    addthis=1; R=[];
    for i=1:size(sStart,1)
        if (Start(k)>=sStart(i))&&(Stop(k)<=sStop(i))
            % // A nested inside S //
            if Lambda(k)>=sLambda(i)
                % -- A is nested inside S of same/lesser Lambda! --
                % retain the smaller cluster, A
                if isempty(find(R==i)), R=[R; i]; end
            else
                addthis=0; break;
            end            
        elseif ((Start(k)<sStart(i))&&(Stop(k)>=sStop(i)))||((Start(k)<=sStart(i))&&(Stop(k)>sStop(i)))
            % // A contains S //
            if Lambda(k)>=sLambda(i)
                % -- A contains S of same/lesser Lambda! --
                addthis=0; break;
            else
                addthis=0; break;
            end
        elseif (length(max(Start(k),sStart(i)):min(Stop(k),sStop(i)))/Nz(k))>0.5
            % // A overlaps significantly with S //
            if Lambda(k)>=sLambda(i)
                % -- A overlaps largely with S of same/lesser Lambda! --                
                % retain the smaller cluster
                if Nz(k)<sNz(i)
                    if isempty(find(R==i)), R=[R; i]; end
                elseif Nz(k)==sNz(i)
                    % if they're both the same size, retain the cluster
                    % with higher Cz
                    if Cz(k)>sCz(i) % sum(x(Start(k):Stop(k)))<sum(x(sStart(i):sStop(i)))
                        if isempty(find(R==i)), R=[R; i]; end
                    else
                        addthis=0; break;
                    end
                else
                    % if A is larger, ignore it
                    addthis=0; break;
                end
            else
                addthis=0; break;
            end
        end
    end
    if addthis==1
        % fprintf('R = %d, length(sLambda) = %d\n',R,length(sLambda));
        sLambda(R)=[]; sStart(R)=[]; sStop(R)=[];  sCz(R)=[]; sNz(R)=[];
        sLambda=[sLambda; Lambda(k)]; sStart=[sStart; Start(k)]; sStop=[sStop; Stop(k)]; sCz=[sCz; Cz(k)]; sNz=[sNz; Nz(k)];
    end
    k=k+1;
end
