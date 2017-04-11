function varargout = Fit_Gui_Chop(varargin)
% FIT_GUI_CHOP MATLAB code for Fit_Gui_Chop.fig
%      FIT_GUI_CHOP, by itself, creates a new FIT_GUI_CHOP or raises the existing
%      singleton*.
%
%      H = FIT_GUI_CHOP returns the handle to a new FIT_GUI_CHOP or the handle to
%      the existing singleton*.
%
%      FIT_GUI_CHOP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIT_GUI_CHOP.M with the given input arguments.
%
%      FIT_GUI_CHOP('Property','Value',...) creates a new FIT_GUI_CHOP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fit_Gui_Chop_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fit_Gui_Chop_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fit_Gui_Chop

% Last Modified by GUIDE v2.5 18-Oct-2016 12:27:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fit_Gui_Chop_OpeningFcn, ...
                   'gui_OutputFcn',  @Fit_Gui_Chop_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end
% --- Executes just before Fit_Gui_Chop is made visible.
function Fit_Gui_Chop_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fit_Gui_Chop (see VARARGIN)
% ===========OPENING SETUP================
if isempty(varargin)
[filename,pathname] = uigetfile('temp_*.mat','Select the temp files');
load ([pathname filename]);
handles.inputfile=[pathname filename];
else
load (varargin{1});%load the temp saved environment.
handles.inputfile=varargin{1};
end
axes(handles.axes1) %select axes1 as the plot unit.
        %Plot the original (Pasted from part 1A:
        figTitle = ['GannetFit Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
        % GABA plot
        % find peak of GABA plot... plot residuals above this...
        gabamin = min(real(GABAData(ii,plotbounds)));
        gabamax = max(real(GABAData(ii,plotbounds)));
        resmax = max(residg);
        residg = residg + gabamin - resmax;
        if strcmp(MRS_struct.p.target,'GABA')
            plot(freq(freqbounds),GaussModel_area(GaussModelParam(ii,:),freq(freqbounds)),'r',...
                freq(plotbounds),real(GABAData(ii,plotbounds)), 'b', ...
                freq(freqbounds),residg,'k');
            set(gca,'XLim',[2.6 3.6]);
            set(gca,'XDir','reverse');
        %end the paste
        %++++++++  labels  +++++++++++
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(3,gabamax/4,MRS_struct.p.target);
        set(hgaba,'horizontalAlignment', 'center');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(GABAData(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(GABAData(ii,upperbound:(upperbound+150))));
        hgabares=text(2.8,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(2.8,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(2.8,tailbottom-gabamax/20,'model','Color',[1 0 0]);
        %++++++++ end of labels +++++++++++         
        end
        
     
   set(handles.text16,'String', num2str(GaussModelParam(1)))
   set(handles.text17,'String', num2str(GaussModelParam(2)))
   set(handles.text18,'String', num2str(GaussModelParam(3)))
   set(handles.text19,'String', num2str(GaussModelParam(4)))
   set(handles.text20,'String', num2str(GaussModelParam(5)))
   %set(handles.edit3,'String', num2str(GaussModelParam(3)))
   %set(handles.slideA,'value',1)
   %set(handles.slideB,'value',1)
   
      %Paste output from code 143 at GannetFit.m
           GABAheight = GaussModelParam(ii,1);
        % FitSTD reports the standard deviation of the residuals / gaba HEIGHT
        MRS_struct.out.GABAFitError(ii)  =  100*std(residg)/GABAheight;
        % This sets GabaArea as the area under the curve.
        MRS_struct.out.GABAArea(ii)=GaussModelParam(ii,1)./sqrt(-GaussModelParam(ii,2))*sqrt(pi);
        sigma = ( 1 / (2 * (abs(GaussModelParam(ii,2)))) ).^(1/2);
        MRS_struct.out.GABAFWHM(ii) =  abs( (2* MRS_struct.p.LarmorFreq) * sigma);
        MRS_struct.out.GABAModelFit(ii,:)=GaussModelParam(ii,:);
        
    %End of past
    
set(handles.text7,'String', num2str(MRS_struct.out.GABAArea(ii)));
set(handles.text8,'String', num2str(MRS_struct.out.GABAFWHM(ii)));
set(handles.text9,'String', num2str(GABAheight));
set(handles.text6,'String', num2str(MRS_struct.out.GABAFitError(ii)));
% ======== END OF OPENING SETUP ===============



% Choose default command line output for Fit_Gui_Chop
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fit_Gui_Chop wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = Fit_Gui_Chop_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function edit_upfreq_input_Callback(hObject, eventdata, handles)
% hObject    handle to edit_upfreq_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_upfreq_input as text
%        str2double(get(hObject,'String')) returns contents of edit_upfreq_input as a double
end

% --- Executes during object creation, after setting all properties.
function edit_upfreq_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_upfreq_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton_v.
function pushbutton_v_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Visulazation button
load (handles.inputfile)
% Get ppm boundary from edits windows
cut_up_freq=str2num(get(handles.edit_upfreq_input, 'String'));
cut_down_freq=str2num(get(handles.edit3, 'String'));
set(handles.text_upfreq_output, 'String', num2str(cut_up_freq));
set(handles.text_downfreq_output, 'String', num2str(cut_down_freq));
Original_Para=GaussModelParam;
% convert ppm to the index
        diff_freq=freq(2)-freq(1);
        cut_low = (cut_up_freq-freq(min(freqbounds)))/diff_freq+1;
        cut_high = (cut_down_freq-freq(min(freqbounds)))/diff_freq+1;
%         cut_low = 350;
%         cut_high = 580;
% chop the index
        freqbounds_chop=[freqbounds(1:cut_low) freqbounds(cut_high:end)];
        freqbounds_orig=freqbounds;
        freqbounds=freqbounds_chop;
        plotbounds=(lowerbound-150):(upperbound+150);
        plotbounds_chop=[plotbounds(1:(cut_low+150)) plotbounds((cut_high+150):end)];
        plotbounds_orig=plotbounds;
        plotbounds=plotbounds_chop;
        
        %Plot
        axes(handles.axes1) %select axes1 as the plot unit.
        %Plot the original (Pasted from part 1A:
        figTitle = ['GannetFit Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
        % GABA plot
        % find peak of GABA plot... plot residuals above this...
        gabamin = min(real(GABAData(ii,plotbounds_orig)));
        gabamax = max(real(GABAData(ii,plotbounds_orig)));
        resmax = max(residg);
        residg = residg + gabamin - resmax;
        if strcmp(MRS_struct.p.target,'GABA')
            plot(freq(freqbounds_orig),GaussModel_area(GaussModelParam(ii,:),freq(freqbounds_orig)),'r',...
                freq(plotbounds_orig),real(GABAData(ii,plotbounds_orig)), 'b', ...
                freq(freqbounds_orig),residg,'k');
            set(gca,'XLim',[2.6 3.6]);
            set(gca,'XDir','reverse');
        %end the paste
        %++++++++  labels  +++++++++++
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(3,gabamax/4,MRS_struct.p.target);
        set(hgaba,'horizontalAlignment', 'center');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(GABAData(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(GABAData(ii,upperbound:(upperbound+150))));
        hgabares=text(2.8,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(2.8,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(2.8,tailbottom-gabamax/20,'model','Color',[1 0 0]);
        %++++++++ end of labels +++++++++++         
        end
        
            YL=ylim;
            try
            l1=line ([cut_down_freq cut_down_freq],[YL(1) real(GABAData(ii,freqbounds_orig(round(cut_high))))],...
                'Color', 'c', 'LineStyle',':','Marker','o','LineWidth',2);
            l2=line ([cut_up_freq cut_up_freq],[YL(1) real(GABAData(ii,freqbounds_orig(round(cut_low))))],...
                'Color', 'c', 'LineStyle',':','Marker','o', 'LineWidth',2);
            catch
                Fit_Warning
            end
end




% --- Executes on button press in pushbutton_c.
function pushbutton_c_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load (handles.inputfile)
% Get ppm boundary from static text windows

cut_up_freq = str2num(get(handles.text_upfreq_output, 'String'));
cut_down_freq = str2num (get(handles.text_downfreq_output, 'String'));
if isempty(cut_up_freq);
    Fit_Warning;
else
Original_Para=GaussModelParam;
% convert ppm to the index
        diff_freq=freq(2)-freq(1);
        cut_low = (cut_up_freq-freq(min(freqbounds)))/diff_freq+1;
        cut_high = (cut_down_freq-freq(min(freqbounds)))/diff_freq+1;
%         cut_low = 350;
%         cut_high = 580;
% chop the index
        freqbounds_chop=[freqbounds(1:cut_low) freqbounds(cut_high:end)];
        freqbounds_orig=freqbounds;
        freqbounds=freqbounds_chop;
        plotbounds=(lowerbound-150):(upperbound+150);
        plotbounds_chop=[plotbounds(1:(cut_low+150)) plotbounds((cut_high+150):end)];
        plotbounds_orig=plotbounds;
        plotbounds=plotbounds_chop;
        % ====Now start the calculation (COPY FROM GANNET FIT)====
        
            [GaussModelParam(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) GaussModel_area(xdummy,ydummy), ...
            GaussModelInit, freq(freqbounds),real(GABAData(ii,freqbounds)), ...
            lb,ub,options);
            residg = -residg;
        if(fit_method == FIT_NLINFIT)
            GaussModelInit = GaussModelParam(ii,:);
            % 1111013 restart the optimisation, to ensure convergence
            for fit_iter = 1:100
                [GaussModelParam(ii,:), residg, J, COVB, MSE] = nlinfit(freq(freqbounds), real(GABAData(ii,freqbounds)), ... % J, COBV, MSE edited in
                    @(xdummy,ydummy) GaussModel_area(xdummy,ydummy), ...
                    GaussModelInit, ...
                    nlinopts);
                MRS_struct.out.fitparams_iter(fit_iter,:,ii) = GaussModelParam(ii,:);
                GaussModelInit = GaussModelParam(ii,:);
                ci = nlparci(GaussModelParam(ii,:), residg,'covar',COVB); %copied over
            end
        end
        % ====END the calculation (COPY FROM GANNET FIT above)====
        
        %Plot
        axes(handles.axes2) %select axes1 as the plot unit.
        %Plot the original (Pasted from part 1A:
        figTitle = ['GannetFit Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
        % GABA plot
        % find peak of GABA plot... plot residuals above this...
        gabamin = min(real(GABAData(ii,plotbounds)));
        gabamax = max(real(GABAData(ii,plotbounds)));
        resmax = max(residg);
        residg = residg + gabamin - resmax;
        if strcmp(MRS_struct.p.target,'GABA')
            plot(freq(freqbounds),GaussModel_area(GaussModelParam(ii,:),freq(freqbounds)),'r',...
            freq(plotbounds),real(GABAData(ii,plotbounds)), 'b', ...
            freq(freqbounds),residg,'k');
%         if strcmp(MRS_struct.p.target,'GABA')
%             plot(freq(freqbounds),GaussModel_area(GaussModelParam(ii,:),freq(freqbounds)),'r',...
%                 freq(plotbounds),real(GABAData(ii,plotbounds)), 'b', ...
%                 freq(freqbounds),residg,'k');
            set(gca,'XLim',[2.6 3.6]);
            set(gca,'XDir','reverse');
        %end the paste
        %++++++++  labels  +++++++++++
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(3,gabamax/4,MRS_struct.p.target);
        set(hgaba,'horizontalAlignment', 'center');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(GABAData(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(GABAData(ii,upperbound:(upperbound+150))));
        hgabares=text(2.8,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(2.8,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(2.8,tailbottom-gabamax/20,'model','Color',[1 0 0]);
        %++++++++ end of labels +++++++++++         
         end
        
            YL=ylim;
            l3=line ([cut_down_freq cut_down_freq],[YL(1) real(GABAData(ii,freqbounds_orig(round(cut_high))))],...
                'Color', 'c', 'LineStyle',':','Marker','o','LineWidth',2);
            l4=line ([cut_up_freq cut_up_freq],[YL(1) real(GABAData(ii,freqbounds_orig(round(cut_low))))],...
                'Color', 'c', 'LineStyle',':','Marker','o', 'LineWidth',2);
            
            
            
    %Outputs
    
   set(handles.text21,'String', num2str(GaussModelParam(1)))
   set(handles.text22,'String', num2str(GaussModelParam(2)))
   set(handles.text23,'String', num2str(GaussModelParam(3)))
   set(handles.text24,'String', num2str(GaussModelParam(4)))
   set(handles.text25,'String', num2str(GaussModelParam(5)))
   %==output=========
        GABAheight = GaussModelParam(ii,1);
        % FitSTD reports the standard deviation of the residuals / gaba HEIGHT
        MRS_struct.out.GABAFitError(ii)  =  100*std(residg)/GABAheight;
        % This sets GabaArea as the area under the curve.
        MRS_struct.out.GABAArea(ii)=GaussModelParam(ii,1)./sqrt(-GaussModelParam(ii,2))*sqrt(pi);
        sigma = ( 1 / (2 * (abs(GaussModelParam(ii,2)))) ).^(1/2);
        MRS_struct.out.GABAFWHM(ii) =  abs( (2* MRS_struct.p.LarmorFreq) * sigma);
        MRS_struct.out.GABAModelFit(ii,:)=GaussModelParam(ii,:);
        
    %End of past
    
set(handles.text13,'String', num2str(MRS_struct.out.GABAArea(ii)));
set(handles.text14,'String', num2str(MRS_struct.out.GABAFWHM(ii)));
set(handles.text15,'String', num2str(GABAheight));
set(handles.text12,'String', num2str(MRS_struct.out.GABAFitError(ii)));

%save the parameters
ifplot=get(handles.radiobutton1,'Value');
if ifplot == 1;
[a1 b1 c1]=fileparts(handles.inputfile);
save ([a1 filesep 'chopped_' b1 c1]);
handles.chopfile=[a1 filesep 'chopped_' b1 c1];
guidata(hObject, handles);
temp2(handles.chopfile);
end


end
end
% --- Executes on button press in pushbutton_r.
function pushbutton_r_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


cla (handles.axes2)

cut_up_freq = str2num(get(handles.text_upfreq_output, 'String'));
cut_down_freq = str2num (get(handles.text_downfreq_output, 'String'));
if ~isempty(cut_up_freq)
set(handles.text30, 'String', num2str(cut_up_freq));
set(handles.text31, 'String', num2str(cut_down_freq));
end
set(handles.edit_upfreq_input, 'String','Pls-input');
set(handles.edit3, 'String','Pls-input');
set(handles.text_upfreq_output, 'String', 'N/A');
set(handles.text_downfreq_output, 'String','N/A');


%Replot==============
% Visulazation button
load (handles.inputfile)
% Get ppm boundary from edits windows
% cut_up_freq=str2num(get(handles.edit_upfreq_input, 'String'));
% cut_down_freq=str2num(get(handles.edit3, 'String'));
% set(handles.text_upfreq_output, 'String', num2str(cut_up_freq));
% set(handles.text_downfreq_output, 'String', num2str(cut_down_freq));
% Original_Para=GaussModelParam;
% convert ppm to the index
%         diff_freq=freq(2)-freq(1);
%         cut_low = (cut_up_freq-freq(min(freqbounds)))/diff_freq+1;
%         cut_high = (cut_down_freq-freq(min(freqbounds)))/diff_freq+1;
% %         cut_low = 350;
% %         cut_high = 580;
% % chop the index
%         freqbounds_chop=[freqbounds(1:cut_low) freqbounds(cut_high:end)];
         freqbounds_orig=freqbounds;
%         freqbounds=freqbounds_chop;
%         plotbounds=(lowerbound-150):(upperbound+150);
%         plotbounds_chop=[plotbounds(1:(cut_low+150)) plotbounds((cut_high+150):end)];
         plotbounds_orig=plotbounds;
%         plotbounds=plotbounds_chop;
%         
        %Plot
        axes(handles.axes1) %select axes1 as the plot unit.
        %Plot the original (Pasted from part 1A:
        figTitle = ['GannetFit Output'];
        set(gcf,'Name',figTitle,'Tag',figTitle, 'NumberTitle','off');
        % GABA plot
        % find peak of GABA plot... plot residuals above this...
        gabamin = min(real(GABAData(ii,plotbounds_orig)));
        gabamax = max(real(GABAData(ii,plotbounds_orig)));
        resmax = max(residg);
        residg = residg + gabamin - resmax;
        if strcmp(MRS_struct.p.target,'GABA')
            plot(freq(freqbounds_orig),GaussModel_area(GaussModelParam(ii,:),freq(freqbounds_orig)),'r',...
                freq(plotbounds_orig),real(GABAData(ii,plotbounds_orig)), 'b', ...
                freq(freqbounds_orig),residg,'k');
            set(gca,'XLim',[2.6 3.6]);
            set(gca,'XDir','reverse');
        %end the paste
        %++++++++  labels  +++++++++++
        %%%%From here on is cosmetic - adding labels (and deciding where to).
        hgaba=text(3,gabamax/4,MRS_struct.p.target);
        set(hgaba,'horizontalAlignment', 'center');
        %determine values of GABA tail (below 2.8 ppm.
        z=abs(MRS_struct.spec.freq-2.79);%2.75
        upperbound=find(min(z)==z);
        tailtop=max(real(GABAData(ii,upperbound:(upperbound+150))));
        tailbottom=min(real(GABAData(ii,upperbound:(upperbound+150))));
        hgabares=text(2.8,min(residg),'residual');
        set(hgabares,'horizontalAlignment', 'left');
        text(2.8,tailtop+gabamax/20,'data','Color',[0 0 1]);
        text(2.8,tailbottom-gabamax/20,'model','Color',[1 0 0]);
        %++++++++ end of labels +++++++++++         
        end
        
%=======end========



end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

end
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(hObject);
delete(handles.figure1);
end