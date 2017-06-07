
function [Lambda,Start,Stop,Cz,Nz]=kscanstat(x,rc_wt_limit,stepsize)

% calculate Kulldorff's scan statistic
N=length(x); 
I=find(x>rc_wt_limit); % center the zones on the rarest codons
LogL=[]; Start=[]; Stop=[]; Cz=[]; Nz=[];
for j=1:length(I)
    % fprintf('\n\nprocessing center position = %d\n',I(j));
    % create zones of increasing "radius" with I(j) as center
    r=1; % start with this radius, then increment by stepsize
    while 1
        % take zone with radius r
        % fprintf('... r=%d\t',r);
        start=max(1,I(j)-r); stop=min(N,I(j)+r); if ((start<1)||(stop>N)), error('start<1 or stop>N'), end
        % fprintf('start=%d,stop=%d\n',start,stop);        
        z=x(start:stop); cz=sum(z); nz=length(z);
        if ((nz/N)>0.5), break, end
        zbar=[];
        if ((start-1)>=1), zbar=[zbar, x(1:start-1)]; end % left side of z
        if ((stop+1)<=N), zbar=[zbar, x(stop+1:N)]; end % right side of z
        if (sum(z)/nz)>(sum(zbar)/(N-nz))
            logL=nz*log(nz/sum(z))+(N-nz)*log((N-nz)/sum(zbar));
            if isnan(logL), error('logL is nan'), end
            % else
            % logL=N*log(N/sum(x));
            % if isnan(logL), error('logL is nan'), end
            LogL=[LogL; logL]; Start=[Start; start]; Stop=[Stop; stop]; Cz=[Cz; cz]; Nz=[Nz; nz];
        end        
        r=r+stepsize;
    end    
end
logL0=N*log(N/sum(x));
Lambda=LogL-logL0;
% sort Lambda in descending order, apply to Start and Stop
[tempLambda,order]=sort(Lambda,'descend'); tempStart=Start(order); tempStop=Stop(order); tempCz=Cz(order); tempNz=Nz(order);

% collect only the distinct clusters
Lambda=[]; Start=[]; Stop=[]; Cz=[]; Nz=[];
i=1; Lambda=[Lambda; tempLambda(i)]; Start=[Start; tempStart(i)]; Stop=[Stop; tempStop(i)]; Cz=[Cz; tempCz(i)]; Nz=[Nz; tempNz(i)]; 
for i=2:length(tempLambda)
    found=0;
    for k=1:length(Lambda)
        if (Lambda(k)==tempLambda(i))&&(Start(k)==tempStart(i))&&(Stop(k)==tempStop(i))&&(Cz(k)==tempCz(i))&&(Nz(k)==tempNz(i))
            found=1; break;
        end
    end
    if (found==0)
        Lambda=[Lambda; tempLambda(i)]; Start=[Start; tempStart(i)]; Stop=[Stop; tempStop(i)]; Cz=[Cz; tempCz(i)]; Nz=[Nz; tempNz(i)]; 
    end
end
