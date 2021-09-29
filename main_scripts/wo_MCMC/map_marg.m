function[estimates,varargout]=map_marg(data,varargin)

  function[stop]=outfunk(estparams,optimValues,state)
   stop=false;
   if strcmp(state,'iter')
      clf
      plot(xk,yk,'k');
      hold on
      scatter(stk,0,[],'gd','filled')
      scatter(estparams,-optimValues.fval,[],'rd','filled')
      xlabel('muM - muF')
      hold off
      drawnow
      disp(estparams)
      disp(optimValues.fval)
      disp(optimValues.iteration)
      disp(optimValues.procedure)
   end
  end
  function[stop]=outfuns(estparams,optimValues,state)
   stop=false;
   if strcmp(state,'iter')
      clf
      plot(xs,ys,'k');
      hold on
      scatter(sts,0,[],'gd','filled')
      scatter(estparams,-optimValues.fval,[],'rd','filled')
      xlabel('sigma')
      hold off
      drawnow
      disp(estparams)
      disp(optimValues.fval)
      disp(optimValues.iteration)
      disp(optimValues.procedure)
   end
  end

mxp=max(data)-min(data);
stkt=0:mxp/50:mxp;
stst=0:mxp/50:mxp;
stk=stkt;
sts=stst;
for i=1:length(stk)
    stk(i)=ppdfk(stkt(i));
end
for i=1:length(sts)
    sts(i)=ppdfs(stst(i));
end
[~,idk]=min(stk);
[~,ids]=min(sts);
stk=stkt(idk);
sts=stst(ids);
% stk=mean(data(data>=mean(data)))-mean(data(data<mean(data)));
% sts=mean([std(data(data>=mean(data))),std(data(data<mean(data)))]);
if ~isempty(varargin)
    xk=varargin{1};
    yk=varargin{2};
    xs=varargin{3};
    ys=varargin{4};
    modelk=@ppdfk;
    models=@ppdfs;
    [estimatesk,valk,exitflagk]=fminsearch(modelk,randn(1).*0.0001+stk,...
        optimset('OutputFcn',@outfunk,'MaxFunEvals',100000,'MaxIter',10000,'TolX',1e-8));
    [estimatess,vals,exitflags]=fminsearch(models,randn(1).*0.0001+sts,...
        optimset('OutputFcn',@outfuns,'MaxFunEvals',100000,'MaxIter',10000,'TolX',1e-8));
else
    modelk=@ppdfk;
    models=@ppdfs;
    exitflagk=0;
    while exitflagk==0
    [estimatesk,valk,exitflagk]=fminsearch(modelk,randn(1).*0.0001+stk,...
        optimset('MaxFunEvals',100000,'MaxIter',300,'TolX',1e-8));
    end
    exitflags=0;
    while exitflags==0
    [estimatess,vals,exitflags]=fminsearch(models,randn(1).*0.0001+sts,...
        optimset('MaxFunEvals',100000,'MaxIter',300,'TolX',1e-8));
    end
end
    function [pdf]=ppdfk(epara)
       pdf=-postk(epara,data);
    end
    function [pdf]=ppdfs(epara)
       pdf=-posts(epara,data);
    end
varargout{1}=-valk;
varargout{2}=exitflagk;
varargout{3}=-vals;
varargout{4}=exitflags;
estimates=[estimatesk,estimatess];
end