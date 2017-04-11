ens=[];
pref=0;   
FILE=cell(0);
PATHN=cell(0);
FILE2=cell(0);
PATHN2=cell(0);

nl=0;
matfile=char('                                      ');
%only want this directory
g=dir;
for i1=1:length(g)
    nametemp=g(i1).name;
    if length(nametemp)>3,
        boolval = 1;
        for i2=1:length(filter_str),
            if isempty(findstr(nametemp,filter_str{i2})),
                boolval = 0;
            end
        end
        if boolval,
            nl=nl+1;
            PATHN{nl}=cd;
            FILE{nl}=nametemp;
            temp = FILE{nl};
            matfile(nl,1:length(temp))=temp;
        end
    end
end

[matfile,I]=sortrows(matfile);
for i1=1:nl,
   FILE2{i1}=FILE{I(i1)};
   PATHN2{i1}=PATHN{I(i1)};
end

FILE=FILE2;
PATHN=PATHN2;
FILE2=cell(0);
PATHN2=cell(0);

fbox=pref*ones([1,length(FILE)]);
n=floor(length(FILE)/30);
taillcol=30*ones([1 n]);
if length(FILE)>n*30
   taillcol=[taillcol length(FILE)-n*30];
end
clear fchoice Abox
fchoice=figure;
set(fchoice,'units','normalized','position',[.02 .06 .96 .86])
numcol=length(taillcol);
for i1=1:numcol
   for j1=1:taillcol(i1)
      numnam=(i1-1)*30+j1;
      Abox(numnam)=uicontrol('style','checkbox','string',deblank(matfile(numnam,:)),'userdata',numnam,...
         'value',fbox(numnam),'units','normalized',...
         'position',[.01+(i1-1)*0.16 .9-(j1-2)*.029 .15 .025],...
         'call','fbox(get(gcbo,''userdata''))=get(gcbo,''value'');');
   end
end
clear matfile
allcall='fbox(:)=0;set(Abox(:),''value'',0)';
revert1box=uicontrol('style','pushbutton','string','all [ ]','call',allcall,...
   'units','normalized','position',[.9 .09 .08 .035]);
allvcall='fbox(:)=1;set(Abox(:),''value'',1)';
revert2box=uicontrol('style','pushbutton','string','all [V]','call',allvcall,...
   'units','normalized','position',[.9 .05 .08 .035]);
strcall='ens2=find(fbox);for k=1:length(ens2),FILE2{k}=FILE{ens2(k)};PATHN2{k}=PATHN{ens2(k)};end;FILE=FILE2;PATHN=PATHN2;close(fchoice);';
okbox=uicontrol('style','pushbutton','string','ok','call',strcall,...
   'units','normalized','position',[.93 .01 .05 .035]);
waitfor(fchoice);