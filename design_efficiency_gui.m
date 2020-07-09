function design_efficiency_gui(hfig,option);
% DESIGN_EFFICIENCY_GUI Graphic user interface to efficiency analyses for fMRI block or event related designs
%
% see also DESIGN_EFFICIENCY
%

% 2008/01
% alfnie@gmail.com

if nargin<1,
	hfig=figure('units','norm','position',[.01,.5,.98,.4],'color','w','name',mfilename,'numbertitle','off','menubar','none','colormap',gray(256));
end
colorb=[1 1 1];
colorc=[.9 .9 .9];%[1,1,.7];
h=get(hfig,'userdata');
if isempty(h) || isempty(h.params),
	if isempty(dir([mfilename,'.mat'])),
		h.params(1)=struct('ncond',3,'conditions',1:3,'blocklengths',[4,4,4],'blockjitters',[0,0,0],'blockdelays',[0,0,0],'tr',2,'hparam',120,'sessionlength',600,'contrast',{{-22,-22}},'sblocklengths',[4,4,4],'dblocklengths',[2,2,2],'iblocklengths',10,'sblockjitters',[0,0,0],'dblockjitters',[1,1,1],'iblockjitters',10,'sblockdelays',[0,0,0],'dblockdelays',[1,1,1],'iblockdelays',10,'iconditionsrepeat',1);
		for n=2:10, h.params(n)=h.params(1);h.params(n).blocklengths=[4,4,4]+2*(n-1); end
		h.n=1;
		h.init=ones(1,length(h.params));
	else,
		load([mfilename,'.mat']);
	end
	h.draw=1;
else,
	h.init=zeros(1,length(h.params));
	h.draw=1;
	switch(option),
		case 'select',
			n=zeros(1,length(h.params));
			for n1=1:length(h.params), n(n1)=get(h.select(n1),'value'); end
			n=find(n); 
			if isempty(n), h.n=1;
			else, n(n==h.n)=[]; if isempty(n), h.n=1; else, h.n=n(1); end; end
		case 'ncond',
			n=str2num(get(h.nconditions,'string'));
            if isempty(n), n=h.params(h.n).ncond; end
			h.params(h.n).conditions=min(n,h.params(h.n).conditions); 
			h.params(h.n).ncond=n;
			h.params(h.n).blocklengths=h.params(h.n).blocklengths(min(length(h.params(h.n).blocklengths),1:h.params(h.n).ncond));
			h.params(h.n).sblocklengths=h.params(h.n).sblocklengths(min(length(h.params(h.n).sblocklengths),1:h.params(h.n).ncond));
			h.params(h.n).dblocklengths=h.params(h.n).dblocklengths(min(length(h.params(h.n).dblocklengths),1:h.params(h.n).ncond));
			h.params(h.n).blockjitters=h.params(h.n).blockjitters(min(length(h.params(h.n).blockjitters),1:h.params(h.n).ncond));
			h.params(h.n).sblockjitters=h.params(h.n).sblockjitters(min(length(h.params(h.n).sblockjitters),1:h.params(h.n).ncond));
			h.params(h.n).dblockjitters=h.params(h.n).dblockjitters(min(length(h.params(h.n).dblockjitters),1:h.params(h.n).ncond));
			h.params(h.n).blockdelays=h.params(h.n).blockdelays(min(length(h.params(h.n).blockdelays),1:h.params(h.n).ncond));
			h.params(h.n).sblockdelays=h.params(h.n).sblockdelays(min(length(h.params(h.n).sblockdelays),1:h.params(h.n).ncond));
			h.params(h.n).dblockdelays=h.params(h.n).dblockdelays(min(length(h.params(h.n).dblockdelays),1:h.params(h.n).ncond));
			for n1=1:2, h.params(h.n).contrast{n1}=min(h.params(h.n).contrast{n1},h.params(h.n).ncond); end
			h.init(h.n)=1;
		case 'conditions',
			n=get(h.conditions,'string');
			if n(1)=='!', h.params(h.n).conditions=str2num(n(2:end));
			else, h.params(h.n).conditions=double(n); end
			h.params(h.n).conditions=h.params(h.n).conditions(:)'-min(h.params(h.n).conditions(:))+1;
			if max(h.params(h.n).conditions)~=h.params(h.n).ncond,
				h.params(h.n).ncond=max(h.params(h.n).conditions);
				set(h.nconditions,'string',num2str(h.params(h.n).ncond));
				set(hfig,'userdata',h);
				design_efficiency_gui(hfig,'ncond');
				return;
			else,
				h.init(h.n)=1;
			end
		case 'blocklengths',
            try, for n1=1:h.params(h.n).ncond, h.params(h.n).blocklengths(n1)=str2num(get(h.blocklength(n1),'string')); end; end
			h.init(h.n)=1;
		case 'blockjitters',
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).blockjitters(n1)=str2num(get(h.blockjitter(n1),'string')); end; end
			h.init(h.n)=1;
		case 'blockdelays',
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).blockdelays(n1)=str2num(get(h.blockdelay(n1),'string')); end; end
			h.init(h.n)=1;
		case 'contrast',
			try, for n1=1:2, h.params(h.n).contrast{n1}=double(regexprep(upper(get(h.contrast(n1),'string')),'[^\*A-Z]',''))-'A'+1; end; end
			h.init(h.n)=1;
		case 'tr',
            n=str2num(get(h.tr,'string'));
			if numel(n)==1, h.params(h.n).tr=n; end
			h.init(h.n)=1;
		case 'hparam',
            n=str2num(get(h.hparam,'string'));
            if numel(n)==1, h.params(h.n).hparam=n; end
			h.init(h.n)=1;
		case 'sessionlength',
            n=str2num(get(h.sessionlength,'string'));
			if numel(n)==1, h.params(h.n).sessionlength=n; end
			h.init(h.n)=1;
		case 'other_conditions',
			h.draw=2;
 		case 'other_blocklengths',
			h.draw=3;
 		case 'other_blockjitters',
			h.draw=4;
 		case 'other_blockdelays',
			h.draw=5;
		case 'other_contrast',
			h.draw=6;
 		case 'done_blocklengths',
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).sblocklengths(n1)=str2num(get(h.sblocklength(n1),'string')); end; end
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).dblocklengths(n1)=str2num(get(h.dblocklength(n1),'string')); end; end
            n=str2num(get(h.iblocklength,'string'));
			if numel(n)==1, h.params(h.n).iblocklengths=n; end
            if h.params(h.n).iblocklengths==1, newn=h.n;
            else newn=length(h.params)+(1:h.params(h.n).iblocklengths);
            end
			for n1=1:numel(newn),
				h.params(newn(n1))=h.params(h.n);
				h.params(newn(n1)).blocklengths=h.params(h.n).sblocklengths+(n1-1)*h.params(h.n).dblocklengths;
				h.init(newn(n1))=1;
			end
 		case 'done_blockjitters',
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).sblockjitters(n1)=str2num(get(h.sblockjitter(n1),'string')); end; end
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).dblockjitters(n1)=str2num(get(h.dblockjitter(n1),'string')); end; end
			n=str2num(get(h.iblockjitter,'string'));
            if numel(n)==1, h.params(h.n).iblockjitters=n; end
            if h.params(h.n).iblockjitters==1, newn=h.n;
            else newn=length(h.params)+(1:h.params(h.n).iblockjitters);
            end
			for n1=1:numel(newn),
				h.params(newn(n1))=h.params(h.n);
				h.params(newn(n1)).blockjitters=h.params(h.n).sblockjitters+(n1-1)*h.params(h.n).dblockjitters;
				h.init(newn(n1))=1;
			end
 		case 'done_blockdelays',
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).sblockdelays(n1)=str2num(get(h.sblockdelay(n1),'string')); end; end
			try, for n1=1:h.params(h.n).ncond, h.params(h.n).dblockdelays(n1)=str2num(get(h.dblockdelay(n1),'string')); end; end
			n=str2num(get(h.iblockdelay,'string'));
            if numel(n)==1, h.params(h.n).iblockdelays=n; end
            if h.params(h.n).iblockdelays==1, newn=h.n;
            else newn=length(h.params)+(1:h.params(h.n).iblockdelays);
            end
			for n1=1:numel(newn),
				h.params(newn(n1))=h.params(h.n);
				h.params(newn(n1)).blockdelays=h.params(h.n).sblockdelays+(n1-1)*h.params(h.n).dblockdelays;
				h.init(newn(n1))=1;
			end
 		case 'done_conditions',
			h.params(h.n).iconditions=get(h.iconditions,'value');
            if isempty(h.params(h.n).iconditions), h.params(h.n).iconditions=1; end
			h.params(h.n).rconditions=get(h.rconditions,'value');
            if isempty(h.params(h.n).rconditions), h.params(h.n).rconditions=0; end
            nrepeat=str2num(get(h.iconditionsrepeat,'string'));
            if isempty(nrepeat), nrepeat=1; end
            if nrepeat==1, newn=h.n;
            else newn=length(h.params)+(1:nrepeat);
            end
			switch(h.params(h.n).iconditions)
				case 1,% fixed order
					h.params(h.n).conditions=[];
					t=0; n=1; while t<=2*h.params(h.n).sessionlength,
						if n==1, 
							new=1;
						else, 
							if h.params(h.n).rconditions, 
								if rem(n,2), new=1+mod(h.params(h.n).conditions(n-2),h.params(h.n).ncond-1);
								else, new=h.params(h.n).ncond; end
							else,
								new=1+mod(n-1,h.params(h.n).ncond);
							end
						end
						h.params(h.n).conditions(n)=new;
						t=t+h.params(h.n).blocklengths(new);
						n=n+1;
					end
				case 2,% random order
                    for n1=1:numel(newn),
                        h.params(newn(n1))=h.params(h.n);
                        
                        h.params(newn(n1)).conditions=[];
                        t=0; n=1; while t<=2*h.params(newn(n1)).sessionlength,
                            if h.params(newn(n1)).rconditions,
                                if rem(n,2), new=ceil((h.params(newn(n1)).ncond-1)*rand);
                                else, new=h.params(newn(n1)).ncond; end
                            else,
                                new=ceil(h.params(newn(n1)).ncond*rand);
                            end
                            h.params(newn(n1)).conditions(n)=new;
                            t=t+h.params(newn(n1)).blocklengths(new);
                            n=n+1;
                        end
                        h.init(newn(n1))=1;
                    end
				case 3,% random order no repetitions
                    for n1=1:numel(newn),
                        h.params(newn(n1))=h.params(h.n);
                        
                        h.params(newn(n1)).conditions=[];
                        t=0; n=1; while t<=2*h.params(newn(n1)).sessionlength,
                            if n==1,
                                if h.params(newn(n1)).rconditions, new=ceil((h.params(newn(n1)).ncond-1)*rand);
                                else, new=ceil((h.params(newn(n1)).ncond)*rand); end
                            else,
                                if h.params(newn(n1)).rconditions,
                                    if rem(n,2), new=1+mod(h.params(newn(n1)).conditions(n-2)+ceil((h.params(newn(n1)).ncond-2)*rand)-1,h.params(newn(n1)).ncond-1);
                                    else, new=h.params(newn(n1)).ncond; end
                                else,
                                    new=1+mod(h.params(newn(n1)).conditions(n-1)+ceil((h.params(newn(n1)).ncond-1)*rand)-1,h.params(newn(n1)).ncond);
                                end
                            end
                            h.params(newn(n1)).conditions(n)=new;
                            t=t+h.params(newn(n1)).blocklengths(new);
                            n=n+1;
                        end
                        h.init(newn(n1))=1;
                    end
			end
			h.init(h.n)=1;
		case 'delete',
			if length(h.params)>1,
				h.params=h.params([1:h.n-1,h.n+1:end]);
				h.E=h.E(:,[1:h.n-1,h.n+1:end]);
				h.Efir=h.Efir(:,[1:h.n-1,h.n+1:end]);
				h.n=min(h.n,length(h.params));
			end
		case 'deleteall',
			if length(h.params)>1,
				h.params=h.params(h.n);
				h.E=h.E(:,h.n);
				h.Efir=h.Efir(:,h.n);
				h.n=1;
			end
		case 'deleteallnewer',
			if length(h.params)>1,
				h.params=h.params(1:h.n);
				h.E=h.E(:,1:h.n);
				h.Efir=h.Efir(:,1:h.n);
				h.n=min(h.n,length(h.params));
			end
		case 'new',
			h.params(end+1)=h.params(h.n);
			h.E(:,end+1)=h.E(:,h.n);
			h.Efir(:,end+1)=h.Efir(:,h.n);
			h.n=length(h.params);
		case 'save',
			save([mfilename,'.mat'],'h');
            fprintf('Results saved to file %s\n',fullfile(pwd,[mfilename,'.mat']));
			return;
		case 'load',
			load([mfilename,'.mat'],'h');
            fprintf('Results loaded from file %s\n',fullfile(pwd,[mfilename,'.mat']));
	end
end

clf;
h.nconditions=uicontrol('units','norm','position',[.025,.85,.1,.08],'style','text','string','Number of conditions','backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
h.nconditions=uicontrol('units','norm','position',[.125,.85,.025,.08],'style','edit','string',num2str(h.params(h.n).ncond),'horizontalalignment','left','tooltipstring','Number of different conditions','callback',[mfilename,'(gcbf,''ncond'')']);
h.conditions=uicontrol('units','norm','position',[.03,.75,.075,.08],'style','text','string','- order','backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
h.conditions=uicontrol('units','norm','position',[.125,.75,.175,.08],'style','edit','string',char('A'+h.params(h.n).conditions-1),'horizontalalignment','left','tooltipstring','Sequential list of conditions presented. Use a standard syntax (e.g. ABABABAB) or a matlab script syntax preceded by the symbol "!" (e.g. !repmat([1,2],[1,4]))','callback',[mfilename,'(gcbf,''conditions'')']);
uicontrol('units','norm','position',[.31,.75,.02,.08],'style','pushbutton','string','...','callback',[mfilename,'(gcbf,''other_conditions'')'],'tooltipstring','Additional options for specifying the block conditions');
h.blocklength=uicontrol('units','norm','position',[.03,.65,.075,.08],'style','text','string',['- duration'],'backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
for n1=1:h.params(h.n).ncond,
  h.blocklength(n1)=uicontrol('units','norm','position',[.125+.025*(n1-1),.65,.025,.08],'style','edit','string',num2str(h.params(h.n).blocklengths(n1)),'horizontalalignment','left','tooltipstring',['Length in seconds of each block/event of condition ',char('A'+n1-1)],'callback',[mfilename,'(gcbf,''blocklengths'')']);
end
uicontrol('units','norm','position',[.135+.025*h.params(h.n).ncond,.65,.02,.08],'style','pushbutton','string','...','callback',[mfilename,'(gcbf,''other_blocklengths'')'],'tooltipstring','Additional options for searching through a range of possible condition lengths');
h.blockdelay=uicontrol('units','norm','position',[.03,.55,.075,.08],'style','text','string',['- onset delay'],'backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
for n1=1:h.params(h.n).ncond,
  h.blockdelay(n1)=uicontrol('units','norm','position',[.125+.025*(n1-1),.55,.025,.08],'style','edit','string',num2str(h.params(h.n).blockdelays(n1)),'horizontalalignment','left','tooltipstring',['onset delay in seconds of each block/event of condition ',char('A'+n1-1)],'callback',[mfilename,'(gcbf,''blockdelays'')']);
end
uicontrol('units','norm','position',[.135+.025*h.params(h.n).ncond,.55,.02,.08],'style','pushbutton','string','...','callback',[mfilename,'(gcbf,''other_blockdelays'')'],'tooltipstring','Additional options for searching through a range of possible condition onset delays');
h.blockjitter=uicontrol('units','norm','position',[.03,.45,.075,.08],'style','text','string',['- onset jitter'],'backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
for n1=1:h.params(h.n).ncond,
  h.blockjitter(n1)=uicontrol('units','norm','position',[.125+.025*(n1-1),.45,.025,.08],'style','edit','string',num2str(h.params(h.n).blockjitters(n1)),'horizontalalignment','left','tooltipstring',['jitter (random onset delay) in seconds of each block/event of condition ',char('A'+n1-1)],'callback',[mfilename,'(gcbf,''blockjitters'')']);
end
uicontrol('units','norm','position',[.135+.025*h.params(h.n).ncond,.45,.02,.08],'style','pushbutton','string','...','callback',[mfilename,'(gcbf,''other_blockjitters'')'],'tooltipstring','Additional options for searching through a range of possible condition jitters');
h.contrast=uicontrol('units','norm','position',[.025,.35,.1,.08],'style','text','string',['Contrast'],'backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
for n1=1:2,
  h.contrast(n1)=uicontrol('units','norm','position',[.125+.025*(n1-1),.35,.025,.08],'style','edit','string',char(h.params(h.n).contrast{n1}-1+'A'),'horizontalalignment','left','tooltipstring',['Specify the conditions pairs (e.g. A-',char('A'+h.params(h.n).ncond-1),') to be included in the contrast of interest (use * to explore all pairwise comparisons)'],'callback',[mfilename,'(gcbf,''contrast'')']);
end
%uicontrol('units','norm','position',[.11+.025*2,.475,.02,.08],'style','pushbutton','string','...','callback',[mfilename,'(gcbf,''other_blocklengths'')'],'tooltipstring','Additional options for testing multiple contrasts simultaneously');
h.tr=uicontrol('units','norm','position',[.025,.25,.1,.08],'style','text','string','RT','backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
h.tr=uicontrol('units','norm','position',[.125,.25,.025,.08],'style','edit','string',num2str(h.params(h.n).tr),'horizontalalignment','left','tooltipstring','Repetition time (in seconds)','callback',[mfilename,'(gcbf,''tr'')']);
h.hparam=uicontrol('units','norm','position',[.025,.15,.1,.08],'style','text','string','High-pass cutoff','backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
h.hparam=uicontrol('units','norm','position',[.125,.15,.05,.08],'style','edit','string',num2str(h.params(h.n).hparam),'horizontalalignment','left','tooltipstring','High-pass filter cutoff (in seconds)','callback',[mfilename,'(gcbf,''hparam'')']);
h.sessionlength=uicontrol('units','norm','position',[.025,.05,.1,.08],'style','text','string','Session length','backgroundcolor',colorb,'fontweight','bold','horizontalalignment','left','fontname','times');
h.sessionlength=uicontrol('units','norm','position',[.125,.05,.05,.08],'style','edit','string',num2str(h.params(h.n).sessionlength),'horizontalalignment','left','tooltipstring','Total session length (in seconds)','callback',[mfilename,'(gcbf,''sessionlength'')']);

h.delete=uicontrol('units','norm','position',[.62,.05,.06,.1],'style','pushbutton','string','Delete','callback',[mfilename,'(gcbf,''delete'')'],'tooltipstring','Deletes currently selected result');
h.deleteall=uicontrol('units','norm','position',[.68,.05,.06,.1],'style','pushbutton','string','Delete all','callback',[mfilename,'(gcbf,''deleteall'')'],'tooltipstring','Deletes all but currently selected result');
h.deleteallnewer=uicontrol('units','norm','position',[.74,.05,.06,.1],'style','pushbutton','string','Delete newer','callback',[mfilename,'(gcbf,''deleteallnewer'')'],'tooltipstring','Deletes all results newer than (to the right of) the currently selected');
h.new=uicontrol('units','norm','position',[.80,.05,.06,.1],'style','pushbutton','string','New','callback',[mfilename,'(gcbf,''new'')'],'tooltipstring','Duplicates the currently selected result');
h.save=uicontrol('units','norm','position',[.87,.05,.06,.1],'style','pushbutton','string','Save','callback',[mfilename,'(gcbf,''save'')'],'tooltipstring',['Saves all results. They will be loaded automatically the next time ',mfilename,' is called from the current directory']);
h.load=uicontrol('units','norm','position',[.93,.05,.06,.1],'style','pushbutton','string','Load','callback',[mfilename,'(gcbf,''load'')'],'tooltipstring',['Loads last-saved results from current directory']);

set(gcf,'pointer','watch');
if nnz(h.init)>10, hmsg=uicontrol('units','norm','position',[.3 .3 .4 .4],'style','text','string','Computing. Please wait...','backgroundcolor',colorb); drawnow; else hmsg=[]; end

for n=find(h.init),
	h.params(n).nscans=h.params(n).sessionlength/h.params(n).tr;
	if h.params(n).contrast{1}<0, IDX1=(1:h.params(n).ncond)'; else, IDX1=h.params(n).contrast{1}; end
	if h.params(n).contrast{2}<0, IDX2=(1:h.params(n).ncond)'; else, IDX2=h.params(n).contrast{2}; end
	h.params(n).contrastvector=zeros(size(IDX1,1)*size(IDX2,1),h.params(n).ncond);
	k=1;
	for n1=1:size(IDX1,1),
		for n2=1:size(IDX2,1),
			if length(IDX1(n1,:))~=length(IDX2(n2,:)) || any(IDX1(n1,:)~=IDX2(n2,:)),
				h.params(n).contrastvector(k,IDX1(n1,:))=1/size(IDX1,2); h.params(n).contrastvector(k,IDX2(n2,:))=-1/size(IDX2,2);
				k=k+1;
			end
		end
	end
	h.params(n).contrastvector=h.params(n).contrastvector(1:k-1,:);
    [temp,h.params(n).dparams]=design_efficiency('conditions',h.params(n).conditions,'blocklengths',h.params(n).blocklengths,'blockjitters',h.params(n).blockjitters,'blockdelays',h.params(n).blockdelays,'tr',h.params(n).tr,'hparam',h.params(n).hparam,'nscans',h.params(n).nscans,'contrast',h.params(n).contrastvector,'model',[3,1]);
	if size(IDX1,1)*size(IDX2,1)>1, h.E(:,n)=[max(temp(:,end));min(temp(:,end))]; h.Efir(:,n)=[max(temp(:,1));min(temp(:,1))]; 
    elseif numel(temp)==2, h.E(:,n)=[temp(:,end);temp(:,end)]; h.Efir(:,n)=[temp(:,1);temp(:,1)]; 
    else h.E(:,n)=[nan;nan]; h.Efir(:,n)=[nan;nan];
    end
end
set(gcf,'pointer','arrow');
if ishandle(hmsg), delete(hmsg); end
switch(h.draw),
	case 1, % design matrix plot
		k1=size(h.params(h.n).dparams.X,1);%;h.params(h.n).nscans*h.params(h.n).dparams.tr/h.params(h.n).dparams.tr0;
		h.design=axes('units','norm','position',[.4,.15,.05,.7]);
		h.designimage=imagesc(h.params(h.n).dparams.X);xlabel('Design matrix');set(gca,'xtick',1:h.params(h.n).ncond,'ylim',[1,k1],'ytick',[1,k1],'yticklabel',{'First scan','Last scan'},'fontsize',6);%set(gca,'xtick',1:h.params(h.n).ncond,'ytick',[1,k1],'yticklabel',{'First scan','Last scan'},'fontsize',6);
		h.line=axes('units','norm','position',[.45,.15,.04,.7]);
		k2=min(k1,ceil(k1/10));%ceil(3*sum(h.params(h.n).blocklengths)/h.params(h.n).dparams.tr0));
		h.lineplot=plot([0,1],[1,1],'k-',[0,1],[(1-k2/k1),.2/.7],'k-'); set(gca,'xlim',[0,1],'ylim',[0,1],'xcolor',colorb,'ycolor',colorb,'xtick',[],'ytick',[],'box','off');
		h.design2=axes('units','norm','position',[.5,.35,.05,.5]);
		h.designimage2=imagesc(1:h.params(h.n).ncond,h.params(h.n).dparams.t((1:k2)),h.params(h.n).dparams.X(1:k2,:));xlabel(strvcat('Design matrix','(close up)'));set(gca,'xtick',1:h.params(h.n).ncond,'fontsize',6);
	case 2, % additional options for generating conditions
		if ~isfield(h.params(h.n),'iconditions'), h.params(h.n).iconditions=1; end
		if ~isfield(h.params(h.n),'rconditions'), h.params(h.n).rconditions=0; end
    	uicontrol('units','norm','position',[.35,.6,.13,.1],'style','text','string',['Multiple results window'],'backgroundcolor','w','fontweight','bold','horizontalalignment','left','fontname','times');
		uicontrol('units','norm','position',[.34,.125,.02+.08+.125,.52],'style','frame','backgroundcolor',colorc);
		uicontrol('units','norm','position',[.35,.37+.205,.08+.125,.05],'style','text','string',['Design types :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		h.iconditions=uicontrol('units','norm','position',[.35,.37,.08+.125,.205],'style','listbox','string',{'Fixed order','Random permutation','Random permutation without repetitions'},'value',h.params(h.n).iconditions);
        uicontrol('units','norm','position',[.41,.20,.13,.05],'style','text','string',['# results :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','right','fontname','times');
        h.iconditionsrepeat=uicontrol('units','norm','position',[.54,.20,.02,.05],'style','edit','string',num2str(h.params(h.n).iconditionsrepeat),'horizontalalignment','left','tooltipstring','number of random design samples (only for ''random permutation'' cases)');
		uicontrol('units','norm','position',[.35,.15,.05,.1],'style','pushbutton','string','Done','callback',[mfilename,'(gcbf,''done_conditions'')'],'tooltipstring',['Imports the selected design in the current result']);
		uicontrol('units','norm','position',[.41,.15,.13,.05],'style','text','string',['interleave rest (',char('A'+h.params(h.n).ncond-1),') block :'],'fontweight','bold','horizontalalignment','right','fontname','times','backgroundcolor',colorc);
		h.rconditions=uicontrol('units','norm','position',[.54,.15,.02,.05],'style','checkbox','value',h.params(h.n).rconditions,'horizontalalignment','center','backgroundcolor',colorc);
	case 3, % additional options for multiple condition lengths
		uicontrol('units','norm','position',[.35,.6,.13,.1],'style','text','string',['Multiple condition lengths window'],'backgroundcolor','w','fontweight','bold','horizontalalignment','left','fontname','times');
		uicontrol('units','norm','position',[.34,.125,.02+.08+.025*h.params(h.n).ncond,.52],'style','frame','backgroundcolor',colorc);
		if ~isfield(h.params(h.n),'sblocklengths'), h.params(h.n).sblocklengths=h.params(h.n).blocklengths; end
		if ~isfield(h.params(h.n),'dblocklengths'), h.params(h.n).dblocklengths=0*h.params(h.n).blocklengths+2; end
		if ~isfield(h.params(h.n),'iblocklengths'), h.params(h.n).iblocklengths=5; end
		h.sblocklength=uicontrol('units','norm','position',[.35,.525,.065,.1],'style','text','string',['starting at :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		for n1=1:h.params(h.n).ncond,
			h.sblocklength(n1)=uicontrol('units','norm','position',[.425+.025*(n1-1),.525,.025,.1],'style','edit','string',num2str(h.params(h.n).sblocklengths(n1)),'horizontalalignment','left','tooltipstring',['Length of each block/event of condition ',char('A'+n1-1),' (in seconds) for the first step']);
		end
		h.dblocklength=uicontrol('units','norm','position',[.35,.4,.065,.1],'style','text','string',['step size :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		for n1=1:h.params(h.n).ncond,
			h.dblocklength(n1)=uicontrol('units','norm','position',[.425+.025*(n1-1),.4,.025,.1],'style','edit','string',num2str(h.params(h.n).dblocklengths(n1)),'horizontalalignment','left','tooltipstring',['Each step will increment the block length of condition ',char('A'+n1-1),' by this amount (in seconds)']);
		end
		uicontrol('units','norm','position',[.35,.275,.065,.1],'style','text','string','# results :','backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		h.iblocklength=uicontrol('units','norm','position',[.425,.275,.025,.1],'style','edit','string',num2str(h.params(h.n).iblocklengths),'horizontalalignment','left','tooltipstring','Number of steps');
		uicontrol('units','norm','position',[.35,.15,.05,.1],'style','pushbutton','string','Done','callback',[mfilename,'(gcbf,''done_blocklengths'')'],'tooltipstring',['Creates a new result in the Design Efficiency window for each step as described above']);
	case 4, % additional options for multiple condition jitters
		uicontrol('units','norm','position',[.35,.6,.13,.1],'style','text','string',['Multiple condition jitters window'],'backgroundcolor','w','fontweight','bold','horizontalalignment','left','fontname','times');
		uicontrol('units','norm','position',[.34,.125,.02+.08+.025*h.params(h.n).ncond,.52],'style','frame','backgroundcolor',colorc);
		if ~isfield(h.params(h.n),'sblockjitters'), h.params(h.n).sblockjitters=0*h.params(h.n).blockjitters; end
		if ~isfield(h.params(h.n),'dblockjitters'), h.params(h.n).dblockjitters=0*h.params(h.n).blockjitters+1; end
		if ~isfield(h.params(h.n),'iblockjitters'), h.params(h.n).iblockjitters=5; end
		h.sblockjitter=uicontrol('units','norm','position',[.35,.525,.065,.1],'style','text','string',['starting at :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		for n1=1:h.params(h.n).ncond,
			h.sblockjitter(n1)=uicontrol('units','norm','position',[.425+.025*(n1-1),.525,.025,.1],'style','edit','string',num2str(h.params(h.n).sblockjitters(n1)),'horizontalalignment','left','tooltipstring',['Length of each block/event of condition ',char('A'+n1-1),' (in seconds) for the first step']);
		end
		h.dblockjitter=uicontrol('units','norm','position',[.35,.4,.065,.1],'style','text','string',['step size :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		for n1=1:h.params(h.n).ncond,
			h.dblockjitter(n1)=uicontrol('units','norm','position',[.425+.025*(n1-1),.4,.025,.1],'style','edit','string',num2str(h.params(h.n).dblockjitters(n1)),'horizontalalignment','left','tooltipstring',['Each step will increment the block length of condition ',char('A'+n1-1),' by this amount (in seconds)']);
		end
		uicontrol('units','norm','position',[.35,.275,.065,.1],'style','text','string','# results :','backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		h.iblockjitter=uicontrol('units','norm','position',[.425,.275,.025,.1],'style','edit','string',num2str(h.params(h.n).iblockjitters),'horizontalalignment','left','tooltipstring','Number of steps');
		uicontrol('units','norm','position',[.35,.15,.05,.1],'style','pushbutton','string','Done','callback',[mfilename,'(gcbf,''done_blockjitters'')'],'tooltipstring',['Creates a new result in the Design Efficiency window for each step as described above']);
    case 5, % additional options for multiple condition delays
		uicontrol('units','norm','position',[.35,.6,.13,.1],'style','text','string',['Multiple condition delays window'],'backgroundcolor','w','fontweight','bold','horizontalalignment','left','fontname','times');
		uicontrol('units','norm','position',[.34,.125,.02+.08+.025*h.params(h.n).ncond,.52],'style','frame','backgroundcolor',colorc);
		if ~isfield(h.params(h.n),'sblockdelays'), h.params(h.n).sblockdelays=0*h.params(h.n).blockdelays; end
		if ~isfield(h.params(h.n),'dblockdelays'), h.params(h.n).dblockdelays=0*h.params(h.n).blockdelays+1; end
		if ~isfield(h.params(h.n),'iblockdelays'), h.params(h.n).iblockdelays=5; end
		h.sblockdelay=uicontrol('units','norm','position',[.35,.525,.065,.1],'style','text','string',['starting at :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		for n1=1:h.params(h.n).ncond,
			h.sblockdelay(n1)=uicontrol('units','norm','position',[.425+.025*(n1-1),.525,.025,.1],'style','edit','string',num2str(h.params(h.n).sblockdelays(n1)),'horizontalalignment','left','tooltipstring',['Length of each block/event of condition ',char('A'+n1-1),' (in seconds) for the first step']);
		end
		h.dblockdelay=uicontrol('units','norm','position',[.35,.4,.065,.1],'style','text','string',['step size :'],'backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		for n1=1:h.params(h.n).ncond,
			h.dblockdelay(n1)=uicontrol('units','norm','position',[.425+.025*(n1-1),.4,.025,.1],'style','edit','string',num2str(h.params(h.n).dblockdelays(n1)),'horizontalalignment','left','tooltipstring',['Each step will increment the block length of condition ',char('A'+n1-1),' by this amount (in seconds)']);
		end
		uicontrol('units','norm','position',[.35,.275,.065,.1],'style','text','string','# results :','backgroundcolor',colorc,'fontweight','bold','horizontalalignment','left','fontname','times');
		h.iblockdelay=uicontrol('units','norm','position',[.425,.275,.025,.1],'style','edit','string',num2str(h.params(h.n).iblockdelays),'horizontalalignment','left','tooltipstring','Number of steps');
		uicontrol('units','norm','position',[.35,.15,.05,.1],'style','pushbutton','string','Done','callback',[mfilename,'(gcbf,''done_blockdelays'')'],'tooltipstring',['Creates a new result in the Design Efficiency window for each step as described above']);
    case 6, % additional options for multiple conditions tested
end
h.efficiency=axes('units','norm','position',[.65,.35,.325,.5]);
maxhe=max(max(h.Efir(:)),max(h.E(:)));
hold(h.efficiency,'on');
n=h.E;
for n1=1:size(n,2), h.efficiencyplot(n1)=patch(n1+.3*[-1,-1,1,1],[n(2,n1),n(1,n1),n(1,n1),n(2,n1)],'k');line(n1+[0,0],[0,n(2,n1)-maxhe*.15],'color',.8*[1,1,1],'linestyle',':');text([n1,n1],[n(1,n1)+maxhe*.1,n(2,n1)-maxhe*.1],{num2str(n(1,n1),'%2.0f'),num2str(n(2,n1),'%2.0f')},'horizontalalignment','center','color','k','fontsize',7); end; 
n=h.Efir;
for n1=1:size(n,2), h.efficiencyplot(n1)=patch(n1+.3*[-1,-1,1,1],[n(2,n1),n(1,n1),n(1,n1),n(2,n1)],'k','facealpha',.5,'facecolor',.8*[1 1 0],'edgecolor',.8*[1 1 0]);line(n1+[0,0],[0,n(2,n1)-maxhe*.15],'color',.8*[1,1,0],'linestyle',':');text([n1,n1],[n(1,n1)+maxhe*.1,n(2,n1)-maxhe*.1],{num2str(n(1,n1),'%2.0f'),num2str(n(2,n1),'%2.0f')},'horizontalalignment','center','color','k','fontsize',7); end; 
ylabel('Design efficiency');set(gca,'xtick',[],'xlim',[0,size(h.E,2)+2],'ylim',[0,max(eps,maxhe*1.25)]);%for n1=1:size(h.E,2),for n2=1:2,if h.E(n2,n1)>0, text(n1,h.E(n2,n1)+max(h.E(:))*(3-2*n2),num2str(h.E(n2,n1),'%2.0f'),'horizontalalignment','center','color','k','fontsize',7); end; end; end; grid on;
hold(h.efficiency,'off');
uicontrol('units','norm','position',[.65 .85 .15 .05],'style','text','string',sprintf('HRF-efficiency %0.1f-%0.1f',h.E(2,h.n),h.E(1,h.n)),'backgroundcolor',colorb,'foregroundcolor','k','fontweight','bold','horizontalalignment','center','fontname','times');
uicontrol('units','norm','position',[.825 .85 .15 .05],'style','text','string',sprintf('FIR-efficiency %0.1f-%0.1f',h.Efir(2,h.n),h.Efir(1,h.n)),'backgroundcolor',colorb,'foregroundcolor',.8*[1 1 0],'fontweight','bold','horizontalalignment','center','fontname','times');

h.select=uicontrol('units','norm','position',[.65,.20,.325,.05],'style','text','string',['Select a result for inspection/editing'],'backgroundcolor',colorb,'fontweight','bold','horizontalalignment','center','fontname','times');
for n1=1:length(h.params),
  h.select(n1)=uicontrol('units','norm','position',[.65+.325/(length(h.params)+2)*n1-.005,.29,.02,.05],'style','checkbox','value',0,'horizontalalignment','center','backgroundcolor','w','callback',[mfilename,'(gcbf,''select'')']);
end
set(h.select(h.n),'value',1);
set(gcf,'userdata',h);
