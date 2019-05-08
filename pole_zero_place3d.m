function varargout = pole_zero_place3d(varargin)
% function pole_zero_place3d(B, A)
% function pole_zero_place3d(B, A, cutoffs)
% function pole_zero_place3d(B, A, gains, fund_frq)
%
% POLE_ZERO_PLACE3D  -- A GUI for filter design using pole-zero placement.
%
% Input Parameters:
%   B       -- feedforward filter coefficients (see help filter) -- optional
%   A       -- feedback filter coefficients (see help filter) -- optional
%         NOTE: If B and A do not define filers with only conjugate pairs,
%               these inputs will be ignored.
%
%     -- For either of the above parameters, you can specify either B or A as
%          empty (i.e., []), and the filter will default to no poles and/or zeros
%
%  * There are three "modes of operation" for this function.  
%    * In the first, we do not specify third or fourth input parameter, and the GUI
%        does not calculate filter evaluation statistics 
%    * In the second, we attempt to design a lowpass filter with high stopband
%        attenuation and low passband ripple.  We use the third input parameter to
%        specify the transition band
%
%   cutoffs -- a two-element vector specifying the starting and ending frequencies
%              of the transition band in Hz
%
%    * In the third, we attempt to design a filter that matches some given frequency
%        response.  We specify the gains at harmonically related frequencies, and
%        a fundamental frequency.
%
%   gains    -- desired transfer function gains at harmonically spaced points
%   fund_frq -- fundamental frequency of the above gains (default: 200 Hz)
%
%
%  >> With this GUI, you can:
%      -- Place pole/zero pairs on the unit circle
%      -- Drag pole/zero pairs around the unit circle
%      -- Move pole/zero pairs in small increments using the arrow keys
%      -- Type in a magnitude and phase for the selected pole/zero pair
%      -- Delete pole/zero pairs  
%      -- Adjust the filter's gain by dragging the frequency response curve
%      -- Display and rotate a 3-D plot of the magnitude of the z-plane
%      -- Switch the frequency response plots from linear units to decibels
%      -- Listen to a harmonic series filtered by your filter (for vowel modeling)
%

% Bad programming practice warning!  This code has not been properly commented at all.
%    Appologies in advance for the insanely illegible spaghetti code...
%
% Some useful tidbits of info:
%   handles.Mode defines several "modes" of operation
%     Mode 0:  Normal operation
%     Mode 1:  Ready to place a zero
%     Mode 2:  Placed/selected zero or pole, waiting for update
%     Mode 4:  Ready to place a pole
%     Mode 5:  Adjusting the gain
%     Mode 6:  Rotating the 3-D plot
%     Mode 8:  Moving poles and zeros using the keyboard
%
% Also note that poles and zeros are treated identically in most places and 
%   are stored in the same arrays.  They can be distinguished using handles.pz
%   0 means zero, 1 means pole.  Similarly, pole/zero pairs and singletons can
%   be distinguished by the handles.pairs variable (pairs - 1, singletons - 0)
%
% The ability to convert pairs to singletons and vice versa was coded into the GUI,
%   but this ability was later deemed confusing.  Thus, now we simply select whether
%   or not to place a singleton or a pair.  The button that allows conversion still
%   exists, but it's hidden on the figure.
%
%  Mark Bartsch, mbartsch@umich.edu  
%  7/9/02

% Written by Mark Bartsch, Winter 2002
% Modification History:
%    8/16/02  --  Added modification history (MB)


if nargin < 1 | ~ischar(varargin{1})  % Create figure
    
    set(0,'ShowHiddenHandles','On')
    h = findobj('Tag','PZ_Place_3D_Figure');
    if ~isempty(h)
        close(h);
    end
    set(0,'ShowHiddenHandles','Off');
    
    fig = openfig(mfilename,'reuse');
        
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    
    handles.fs = 8192;
    default_frq = 200/handles.fs*2*pi;
    default_num_harmonics = floor(pi/default_frq);
    
    
    if nargin < 1 | isempty(varargin{1})
        handles.B = [1];
    else
        handles.B = varargin{1}(:)';
    end
    
    if nargin < 2 | isempty(varargin{2})
        handles.A = [1];
    else
        handles.A = varargin{2}(:)'; 
    end
    [handles.B,handles.A] = eqtflength(handles.B,handles.A);
    
    if ~isreal(handles.B) | ~isreal(handles.A)
        warning('Filter coefficients must be real.');
        handles.B = 1;
        handles.A = 1;
    end
    
    if handles.A(1) == 0
        warning('Leading feedback coefficient must be nonzero.');
        handles.A = 1;
        handles.B = 1;
    end
    
    if nargin < 4
        handles.match_frq = default_frq;
    else
        handles.match_frq = varargin{4};
    end
    
    
    if nargin < 3
        handles.filter_eval = 0;
        set(handles.lin_err_txt,'Visible','Off');
        set(handles.db_err_txt,'Visible','Off');
        set(handles.LinError,'Visible','Off');
        set(handles.dBError,'Visible','Off');
        set(handles.text13,'Visible','Off');
        set(handles.frame3,'Visible','Off');
        
    elseif nargin < 4
        handles.filter_eval = 1;
        handles.cutoff = varargin{3}(1:2);
        
        set(handles.lin_err_txt,'String','Passband Ripple (dB)');
        set(handles.db_err_txt,'String','Stopband Attenuation (dB)');
        
    else
        handles.filter_eval = 2;
        handles.match_vals = varargin{3};
        handles.match_vals((1:length(handles.match_vals))*handles.match_frq > pi) = [];
        
        if size(handles.match_vals,1) > size(handles.match_vals,2)
            handles.match_vals = handles.match_vals';
        end
        
    end
    
    % Update handles structure
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
    set(fig,'HandleVisibility','On');
    
    % Perform a some initialization
    init(fig,[],handles);
    
    set(fig,'HandleVisibility','Callback');
    
    
else % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end
    
end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function init(h, eventdata, handles, varargin)


handles.ax_ext = 1.2;
draw_axes(h,[],handles);

% Store Data in handles

set(gcf,'units','normalized');
set(gcf,'doublebuffer','on');

[z,p,k] = tf2zp(handles.B,handles.A);
z(imag(z) < 0) = [];
p(imag(p) < 0) = [];

handles.zero_list = [z; p];
handles.pz = [zeros(size(z)); ones(size(p))];
handles.pairs = imag(handles.zero_list) ~= 0;

f = find(abs(handles.zero_list) == 0 & handles.pairs == 0 & handles.pz == 0);
while length(f) > 1
    handles.zero_list(f(1)) = [];
    handles.pz(f(1)) = [];
    handles.pairs(f(2)) = 1;
    handles.pairs(f(1)) = [];
    f = find(abs(handles.zero_list) == 0 & handles.pairs == 0 & handles.pz == 0);
end

f = find(abs(handles.zero_list) == 0 & handles.pairs == 0 & handles.pz == 1);
while length(f) > 1
    handles.zero_list(f(1)) = [];
    handles.pz(f(1)) = [];
    handles.pairs(f(2)) = 1;
    handles.pairs(f(1)) = [];
    f = find(abs(handles.zero_list) == 0 & handles.pairs == 0 & handles.pz == 1);
end

handles.zero_handle_list = [];
axes(handles.ZeroPlot);

for i = 1:size(handles.zero_list,1)
    
    if handles.pairs(i)
        hdl = line(real(handles.zero_list(i))*[1 1],imag(handles.zero_list(i))*[1 -1]);
    else
        hdl = line(real(handles.zero_list(i))*[1 1],[0 0]);
    end
    
    set(hdl,'LineStyle','none','LineWidth',1,'MarkerSize',8,'Color',[0 1 0]);
    set(hdl,'ButtonDownFcn',sprintf('pole_zero_place3d(''zero_select'',gcbo,[],guidata(gcbo));'));
%    set(hdl,'EraseMode','xor');
    handles.zero_handle_list = [handles.zero_handle_list hdl];
    
    if handles.pz(i,1) == 0
        set(hdl,'Marker','o');
    else
        set(hdl,'Marker','x');
    end
    
end

axis square
handles.warning_text = text(0,0.2,'WARNING: Unstable Filter!',...
    'FontSize',14,'Color',[1 0 0],'HorizontalAlignment','center','Visible','off');

handles.multiplicity = text(.1,.1,'0','FontSize',8,'HorizontalAlignment','center','Visible','off');

set(gcf,'doublebuffer','on');
set(handles.Gain,'string',num2str(k));
handles.selected_zero = 0;
handles.Mode = 0;


axes(handles.MagResponse);

[fr,w] = freqz(handles.B,handles.A,1024);
%w = w/2/pi*handles.fs;

if handles.filter_eval == 2
    if length(handles.match_vals) < 50
        handles.match_plot = stem(handles.match_frq*(1:length(handles.match_vals)),handles.match_vals,'r');
    else
        handles.match_plot = plot(handles.match_frq*(1:length(handles.match_vals)),handles.match_vals,'r');
    end
elseif handles.filter_eval == 1
    handles.passband = find(w < handles.cutoff(1));
    handles.stopband = find(w > handles.cutoff(2));
    
    handles.pass_min = min(abs(fr(handles.passband)));
    handles.pass_max = max(abs(fr(handles.passband)));
    handles.stop_max = max(abs(fr(handles.stopband)));
    
    ax = axis;
    cutoff = handles.cutoff;
    handles.pmax_line = line([0 cutoff(1) cutoff(1)],[handles.pass_max handles.pass_max ax(4)],'Color',[1 0 0]); 
    handles.pmin_line = line([0 cutoff(1) cutoff(1)],[handles.pass_min handles.pass_min 0],'Color',[1 0 0]); 
    handles.smin_line = line([pi cutoff(2) cutoff(2)],[handles.stop_max handles.stop_max ax(4)],'Color',[1 0 0]); 
    
    if get(handles.LinearButton,'value')
        set(handles.LinError,'string',num2str(abs(20*log10(handles.pass_max/handles.pass_min))));
        set(handles.dBError,'string',num2str(20*log10(handles.pass_max/handles.stop_max)));
    else
        set(handles.LinError,'string',num2str(abs(handles.pass_max - handles.pass_min)));
        set(handles.dBError,'string',num2str(handles.pass_max - handles.stop_max));            
    end
    
end

hold on;
handles.fr = plot(w,abs(fr),'b','EraseMode','xor', ...
    'ButtonDownFcn',sprintf('pole_zero_place3d(''change_gain'',gcbo,[],guidata(gcbo));'));
hold off;
y_lab_h = ylabel('|H(\omega)|');
xlabel('\omega (Discrete radian frequency)');


axes(handles.surf_plot)

handles.num_pts = 46;
lspace = linspace(-handles.ax_ext,handles.ax_ext,handles.num_pts);
[x,y] = meshgrid(lspace,lspace);
handles.z = x+j*y;
handles.surf = mesh(x,y,ones(size(x)));

set(handles.surf_plot,'ButtonDownFcn',...
    'pole_zero_place3d(''surf_plot_ButtonDownFcn'',gcbo,[],guidata(gcbo));');
set(handles.surf,'ButtonDownFcn',...
    'pole_zero_place3d(''surf_plot_ButtonDownFcn'',gcbo,[],guidata(gcbo));');

handles.ce = exp(-j*2*pi*(0:.001:1));

handles.ce1 = line(real(handles.ce),imag(handles.ce),1.01*ones(size(handles.ce)));
set(handles.ce1,'Color',[1 0 0]);
set(handles.ce1,'LineWidth',2);

handles.ce2 = line(real(handles.ce),imag(handles.ce),0.99*ones(size(handles.ce)));
set(handles.ce2,'Color',[1 0 0]);
set(handles.ce2,'LineWidth',2);

xlabel('Real(z)');
ylabel('Imag(z)');
zlabel('|H(z)|');

handles.upplus = 2;
handles.log_upplus = 30;

if ~get(handles.threeDon,'Value')
    set(handles.surf,'Visible','Off');
    set(handles.ce1,'Visible','Off');
    set(handles.ce2,'Visible','Off');
end

guidata(h,handles);

slow_update_frq_resp(h,[],handles);


% --------------------------------------------------------------------
function slow_update_frq_resp(h, eventdata, handles, varargin);


origin_zeros = find(abs(handles.zero_list)==0 & handles.pz == 0);
origin_poles = find(abs(handles.zero_list)==0 & handles.pz == 1);

num_zeros = length(origin_zeros) + sum(handles.pairs(origin_zeros));
num_poles = length(origin_poles) + sum(handles.pairs(origin_poles));
pz = num_zeros - num_poles;

if (num_poles > 0 & num_zeros > 0)

    pairs = find(abs(handles.zero_list)==0 & handles.pairs == 1);
    singles = find(abs(handles.zero_list)==0 & handles.pairs == 0);
    if pz < 0
        num_pairs = floor(abs(pz)/2);
        handles.pz(pairs(1:num_pairs)) = 1;
        set(handles.zero_handle_list(pairs(1:num_pairs)),'Marker','x');
        pairs(1:num_pairs) = [];

        num_singles = abs(pz) - 2*num_pairs;
        handles.pz(singles(1:num_singles)) = 1;
        set(handles.zero_handle_list(singles(1:num_singles)),'Marker','x');
        singles(1:num_singles) = [];
    else
        num_pairs = floor(abs(pz)/2);
        handles.pz(pairs(1:num_pairs)) = 0;
        set(handles.zero_handle_list(pairs(1:num_pairs)),'Marker','o');
        pairs(1:num_pairs) = [];

        num_singles = abs(pz) - 2*num_pairs;
        handles.pz(singles(1:num_singles)) = 0;
        set(handles.zero_handle_list(singles(1:num_singles)),'Marker','o');
        singles(1:num_singles) = [];
    end
    
    delete(handles.zero_handle_list([pairs; singles]));
    handles.zero_handle_list([pairs; singles]) = [];
    handles.zero_list([pairs; singles]) = [];
    handles.pz([pairs; singles]) = [];
    handles.pairs([pairs; singles]) = [];
end

num_origin = sum((abs(handles.zero_list)== 0)) + sum(handles.pairs(abs(handles.zero_list) == 0));
if num_origin > 0
    set(handles.multiplicity,'String',num2str(num_origin),'Visible','On');
else
    set(handles.multiplicity,'Visible','Off');
end        


axes(handles.MagResponse);

num_z = sum((handles.pz == 0).*(1+handles.pairs));
num_orig_z = sum((handles.pz == 0).*(abs(handles.zero_list) == 0).*(1+handles.pairs)); 

num_p = sum((handles.pz == 1).*(1+handles.pairs));
num_orig_p = sum((handles.pz == 1).*(abs(handles.zero_list) == 0).*(1+handles.pairs)); 

set(handles.NumZeros,'string',num2str(num_z - num_orig_z));
set(handles.NumOriginZeros,'String',num2str(num_orig_z));
set(handles.NumPoles,'string',num2str(num_p - num_orig_p));
set(handles.NumOriginPoles,'string',num2str(num_orig_p));

dif = num_z - num_p;
if dif < 0
    set(handles.NumInfZeros,'string',num2str(abs(dif)));
    set(handles.NumInfPoles,'string','0');
elseif dif > 0
    set(handles.NumInfPoles,'string',num2str(abs(dif)));
    set(handles.NumInfZeros,'string','0');
end
    
zlist = handles.zero_list(find(handles.pz == 0));
zlist = zlist - j*imag(zlist).*(1-handles.pairs(find(handles.pz == 0)));
zlist = [zlist; conj(handles.zero_list(find(handles.pz == 0 & handles.pairs)))];
plist = handles.zero_list(find(handles.pz == 1));
plist = plist - j*imag(plist).*(1-handles.pairs(find(handles.pz == 1)));
plist = [plist; conj(handles.zero_list(find(handles.pz == 1 & handles.pairs)))];

if any((abs(plist) >= 1))
    set(handles.warning_text,'Visible','on');
else
    set(handles.warning_text,'Visible','off');
end

[handles.B,handles.A] = zp2tf(zlist,plist,str2double(get(handles.Gain,'string')));

[fr,w] = freqz(handles.B,handles.A,1024);

fr = abs(fr);
if get(handles.LinearButton,'value')
    set(handles.fr,'YData',abs(fr));
else
    warning off
    fr = 20*log10(fr);
    %    warning on
    %    m = mean(fr);
    %    fr(fr < m-50) = m-50;
    %    fr(fr > m+50) = m+50;
    set(handles.fr,'YData',fr);
end

if handles.filter_eval == 1
    set(handles.pmax_line,'YData',max(fr)*[1 1 1]);
    set(handles.pmin_line,'YData',max(fr)*[1 1 1]);
    set(handles.smin_line,'YData',max(fr)*[1 1 1]);
end

axis tight;
axis autoy;

if handles.filter_eval == 2
    
    if length(handles.match_vals) == 1
        H = freqz(handles.B,handles.A,[handles.match_frq 1]);
        H = H(1);
    else
        H = freqz(handles.B,handles.A,(handles.match_frq*(1:length(handles.match_vals))));
    end
    lin_err = sqrt(mean((handles.match_vals - abs(H)).^2));
    warning off
    tmp = 20*log10(handles.match_vals) - 20*log10(abs(H));
    tmp(isinf(tmp) | isnan(tmp)) = [];
    db_err = sqrt(mean(tmp.^2));
    warning on
    
    set(handles.LinError,'string',num2str(lin_err));
    set(handles.dBError,'string',num2str(db_err));
    
elseif handles.filter_eval == 1
    
    handles.passband = find(w < handles.cutoff(1));
    handles.stopband = find(w > handles.cutoff(2));
    
    handles.pass_min = min(fr(handles.passband));
    handles.pass_max = max(fr(handles.passband));
    handles.stop_max = max(fr(handles.stopband));
    
    ax = axis;
    set(handles.pmax_line,'YData',[handles.pass_max handles.pass_max ax(4)]);
    set(handles.pmin_line,'YData',[handles.pass_min handles.pass_min ax(3)]);
    set(handles.smin_line,'YData',[handles.stop_max handles.stop_max ax(4)]);
    
    if get(handles.LinearButton,'value')
        set(handles.LinError,'string',num2str(abs(20*log10(handles.pass_max/handles.pass_min))));
        set(handles.dBError,'string',num2str(20*log10(handles.pass_max/handles.stop_max)));
    else
        set(handles.LinError,'string',num2str(abs(handles.pass_max - handles.pass_min)));
        set(handles.dBError,'string',num2str(handles.pass_max - handles.stop_max));            
    end
    
end

if get(handles.threeDon,'Value')
    
    axes(handles.surf_plot)
    cevals = abs(polyval(handles.B,handles.ce)./polyval(handles.A,handles.ce));
    vals = abs(polyval(handles.B(end:-1:1),1./handles.z)./polyval(handles.A(end:-1:1),1./handles.z));
    
    if get(handles.LinearButton,'value')
        vals(vals > max(cevals)+handles.upplus) = max(cevals)+handles.upplus;
        zmax = max([max(max(vals)),cevals+0.01]);
        zmin = 0;
    else
        warning off
        cevals = 20*log10(cevals);
        warning on
        m = mean(cevals(~isinf(cevals) & ~isnan(cevals)));
        cevals(cevals > 50+m ) = m+50;
        cevals(cevals < m-50) = m-50;
        vals = 20*log10(abs(vals));
        vals(vals > max(cevals)+handles.log_upplus) = max(cevals)+handles.log_upplus;
        vals(vals < min(cevals)-handles.log_upplus) = min(cevals)-handles.log_upplus;
        zmax = max([max(max(vals)),cevals+1]);
        zmin = min([min(min(vals)),cevals-1]);
    end
    
    c = campos;
    set(handles.surf,'ZData',vals);
    set(handles.surf,'CData',vals);
    set(handles.ce1,'ZData',cevals+0.01);
    set(handles.ce2,'ZData',cevals-0.01);
    
    axis([-handles.ax_ext,handles.ax_ext,-handles.ax_ext,handles.ax_ext,zmin,zmax]);
    
    if (c(1) ~= .5 & c(2) ~=.5)
        c2 = campos;
        campos([c(1:2) c2(3)]);
    end
    
end

guidata(h,handles);

% --------------------------------------------------------------------

function draw_axes(h, eventdata, handles, varargin)

axes(handles.ZeroPlot);

ce = exp(j*2*pi*(0:.01:1));

h = line(real(ce),imag(ce));
set(h,'linestyle',':','color',[0 0 0],'HitTest','off')

h = line([0 0],[-handles.ax_ext handles.ax_ext]);
set(h,'linestyle',':','color',[0 0 0],'HitTest','off')
h = line([-handles.ax_ext handles.ax_ext],[0 0]);
set(h,'linestyle',':','color',[0 0 0],'HitTest','off')
axis([-handles.ax_ext handles.ax_ext -handles.ax_ext handles.ax_ext]);

xlabel('Real(z)');
ylabel('Imag(z)');


% --------------------------------------------------------------------
function varargout = ZeroPlot_CreateFcn(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = ZeroMag_Callback(h, eventdata, handles, varargin)

if handles.selected_zero > 0
    try
        newabs = eval(get(h,'string'));
        zeros = handles.zero_list(handles.selected_zero);
        if handles.pairs(handles.selected_zero)
            zeros = zeros./abs(zeros)*newabs;
            if imag(zeros) < 0
                zeros = conj(zeros);
            end
            handles.zero_list(handles.selected_zero) = zeros;
            set(handles.zero_handle_list(handles.selected_zero),'XData',real(zeros)*[1 1], 'YData', imag(zeros)*[1 -1]);
        else
            zeros = newabs*sign(real(zeros)) + j*abs(tan(angle(zeros))*newabs);
            handles.zero_list(handles.selected_zero) = zeros;
            set(handles.zero_handle_list(handles.selected_zero),'XData',real(zeros), 'YData', 0);
            zeros = real(zeros);
        end
            
        set(h,'string',num2str(abs(newabs)));
        set(handles.ZeroAngle,'string',num2str(angle(zeros)));        
        slow_update_frq_resp(h,[],handles);
    catch
        set(h,'string',num2str(abs(handles.zero_list(handles.selected_zero,1))));
    end
else
    set(h,'string','');
end

% --------------------------------------------------------------------
function varargout = ZeroAngle_Callback(h, eventdata, handles, varargin)

if handles.selected_zero > 0
    try
        newang = eval(get(h,'string'));
        newang = mod(newang,2*pi);
        if newang > pi
            newang = newang - 2*pi;
        end
        set(h,'string',num2str(newang));
        zeros = handles.zero_list(handles.selected_zero);
        if handles.pairs(handles.selected_zero)
            zeros = abs(zeros).*exp(j*[newang]);
            if imag(zeros) < 0
                zeros = conj(zeros);
            end
            set(handles.zero_handle_list(handles.selected_zero),'XData',real(zeros)*[1 1], 'YData', imag(zeros)*[1 -1]);
        else
            zeros = real(zeros)*(1 + j*abs(tan(newang)));
            set(handles.zero_handle_list(handles.selected_zero),'XData',real(zeros), 'YData', 0);
            set(h,'string',num2str(angle(real(zeros))));
        end
        handles.zero_list(handles.selected_zero) = zeros;
        
        slow_update_frq_resp(h,[],handles);
        
    catch
        set(h,'string',num2str(abs(handles.zero_list(handles.selected_zero))));
    end
else
    set(h,'string','');
end

% --------------------------------------------------------------------
function varargout = DeleteZero_Callback(h, eventdata, handles, varargin)

if handles.selected_zero > 0
    handles.zero_list(handles.selected_zero) = 0;
    
    set(handles.zero_handle_list(handles.selected_zero),'XData',[0 0]);
    set(handles.zero_handle_list(handles.selected_zero),'YData',[0 -0]);
    set(handles.zero_handle_list(handles.selected_zero),'LineWidth',1);
    
    handles.selected_zero = 0;
    
    slow_update_frq_resp(h,[],handles);
    
end


% --------------------------------------------------------------------
function varargout = MagResponse_CreateFcn(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = AddZero_Callback(h, eventdata, handles, varargin)

if eventdata == 0
    if get(handles.AddZero,'value')
        handles.Mode = 1;
        set(handles.AddPole,'value',0)
    else
        handles.Mode = 0;
    end
else
    if get(handles.AddPole,'value')
        handles.Mode = 4;
        set(handles.AddZero,'value',0)
    else
        handles.Mode = 0;
    end
end

guidata(h,handles);


% --------------------------------------------------------------------
function varargout = ZeroPlot_ButtonDownFcn(h, eventdata, handles, varargin)

if handles.Mode == 1 | handles.Mode == 4
    pt = get(handles.ZeroPlot,'currentpoint');
    pt = pt(1,1:2);
    pt(2) = abs(pt(2));
    handles.zero_list = [pt(1) + j*pt(2); handles.zero_list; 0];
    handles.pairs = [1; handles.pairs; 1];
    handles.selected_zero = 1;
    set(handles.pair_cb,'value',1);
    
    axes(handles.ZeroPlot);
    hdl = line([pt(1) pt(1)],[pt(2) -pt(2)]);
    set(hdl,'LineStyle','none','LineWidth',2,'MarkerSize',8,'Color',[0 1 0]);
    set(hdl,'ButtonDownFcn',sprintf('pole_zero_place3d(''zero_select'',gcbo,[],guidata(gcbo));'));
    for ind = 1:length(handles.zero_handle_list)
        set(handles.zero_handle_list(ind),'LineWidth',1);
    end
    
    hdl2 = line([0 0],[0 -0]);
    set(hdl2,'LineStyle','none','LineWidth',1,'MarkerSize',8,'Color',[0 1 0]);
    set(hdl2,'ButtonDownFcn',sprintf('pole_zero_place3d(''zero_select'',gcbo,[],guidata(gcbo));'));
    
    if handles.Mode == 1
        set(hdl,'Marker','o');
        set(hdl2,'Marker','x');
        handles.pz = [0; handles.pz; 1];
    else
        set(hdl,'Marker','x');
        set(hdl2,'Marker','o');
        handles.pz = [1; handles.pz; 0];
    end
    
    handles.Mode = 2;
    set(handles.AddZero,'value',0);
    set(handles.AddPole,'value',0);
    
    handles.zero_handle_list = [hdl handles.zero_handle_list hdl2]; 

    zero = handles.zero_list(handles.selected_zero);
    pt = get(handles.ZeroPlot,'currentpoint');
    if pt(1,2) < 0
        zero = conj(zero);
    end

    if handles.pairs(handles.selected_zero)
        set(handles.ZeroMag,'string',num2str(abs(zero)));
        set(handles.ZeroAngle,'string',num2str(angle(zero)));   
    else
        set(handles.ZeroMag,'string',num2str(abs(real(zero))));
        set(handles.ZeroAngle,'string',num2str(angle(real(zero))));   
    end        
end

slow_update_frq_resp(h,[],handles);
  


% --------------------------------------------------------------------
function varargout = zero_select(h, eventdata, handles, varargin)

if handles.Mode == 0
    
    for ind = 1:length(handles.zero_handle_list)
        set(handles.zero_handle_list(ind),'LineWidth',1);
    end
    
    handles.selected_zero = find(handles.zero_handle_list == h);
    set(h,'LineWidth',2);
    
    handles.Mode = 2;
    
    zero = handles.zero_list(handles.selected_zero);
    pt = get(handles.ZeroPlot,'currentpoint');
    if pt(1,2) < 0
        zero = conj(zero);
    end

    if handles.pairs(handles.selected_zero)
        set(handles.ZeroMag,'string',num2str(abs(zero)));
        set(handles.ZeroAngle,'string',num2str(angle(zero)));   
    else
        set(handles.ZeroMag,'string',num2str(abs(real(zero))));
        set(handles.ZeroAngle,'string',num2str(angle(real(zero))));   
    end        

    set(handles.pair_cb,'Value',handles.pairs(handles.selected_zero));
    
end

guidata(h,handles);


% --------------------------------------------------------------------
function varargout = change_gain(h, eventdata, handles, varargin)

if handles.Mode == 0
    
    handles.Mode = 5;
    
    pt = get(handles.MagResponse,'currentpoint');
    handles.init_point = pt(3);
    handles.init_gain = str2double(get(handles.Gain,'string'));
    
end

guidata(h,handles);


% --------------------------------------------------------------------
function varargout = figure1_WindowButtonUpFcn(h, eventdata, handles, varargin)

if handles.Mode == 2
    handles.Mode = 0;

    zero = handles.zero_list(handles.selected_zero);
    pt = get(handles.ZeroPlot,'currentpoint');
    if pt(1,2) < 0
        zero = conj(zero);
    end
    
    if handles.pairs(handles.selected_zero)
        set(handles.ZeroMag,'string',num2str(abs(zero)));
        set(handles.ZeroAngle,'string',num2str(angle(zero)));   
    else
        set(handles.ZeroMag,'string',num2str(abs(real(zero))));
        set(handles.ZeroAngle,'string',num2str(angle(real(zero))));   
    end        
    
    slow_update_frq_resp(h,[],handles);
    
elseif handles.Mode == 5
    
    handles.Mode = 0;
    slow_update_frq_resp(h,[],handles);
    
elseif handles.Mode == 6
    
    handles.Mode = 0;
    guidata(h,handles);
    
end



% --------------------------------------------------------------------
function varargout = Gain_Callback(h, eventdata, handles, varargin)

if isnan(str2double(get(h,'string')))
    set(h,'string','1.0');
end

handles.B = zp2tf(handles.zero_list(:),1,str2double(get(handles.Gain,'string')));
slow_update_frq_resp(h,[],handles);


% --------------------------------------------------------------------
function varargout = figure1_WindowButtonMotionFcn(h, eventdata, handles, varargin)

if handles.Mode == 2 | handles.Mode == 5 | handles.Mode == 8
    
    if handles.Mode == 2 | handles.Mode == 8
        if handles.Mode == 2
            pt = get(handles.ZeroPlot,'currentpoint');
            
        else    % handles.Mode == 8;
            pt = handles.key_pt;
            handles.Mode = 0;    
        end
        
        pt = pt(1,1:2);
        pt(2) = abs(pt(2));
        
        if handles.pairs(handles.selected_zero)
            set(handles.zero_handle_list(handles.selected_zero),'XData',[pt(1) pt(1)]);
            set(handles.zero_handle_list(handles.selected_zero),'YData',[pt(2) -pt(2)]);
            handles.zero_list(handles.selected_zero) = [pt(1) + j*pt(2)];
        else
            set(handles.zero_handle_list(handles.selected_zero),'XData',pt(1),'YData',0);
            handles.zero_list(handles.selected_zero) = pt(1)*(1 + j*abs(tan(angle(handles.zero_list(handles.selected_zero)))));
        end            

        zero = handles.zero_list(handles.selected_zero);
        if pt(1,2) < 0
            zero = conj(zero);
        end
        
        if handles.pairs(handles.selected_zero)
            set(handles.ZeroMag,'string',num2str(abs(zero)));
            set(handles.ZeroAngle,'string',num2str(angle(zero)));   
        else
            set(handles.ZeroMag,'string',num2str(abs(real(zero))));
            set(handles.ZeroAngle,'string',num2str(angle(real(zero))));   
        end        

        
        num_origin = sum((abs(handles.zero_list)== 0)) + sum(handles.pairs(abs(handles.zero_list) == 0));
        if num_origin > 0
            set(handles.multiplicity,'String',num2str(num_origin),'Visible','On');
        else
            set(handles.multiplicity,'Visible','Off');
        end        
       
    else
        
        pt = get(handles.MagResponse,'currentpoint');
        if get(handles.LinearButton,'value')
            ratio = pt(3)/handles.init_point;
            new_gain = abs(handles.init_gain*ratio);
        else
            dif = pt(3) - handles.init_point;
            ratio = 10.^(dif/20);
            new_gain = abs(handles.init_gain*ratio);
        end        
        
        set(handles.Gain,'string',num2str(new_gain));
        
    end
    
    zlist = handles.zero_list(find(handles.pz == 0));
    zlist = zlist - j*imag(zlist).*(1-handles.pairs(find(handles.pz == 0)));
    zlist = [zlist; conj(handles.zero_list(find(handles.pz == 0 & handles.pairs)))];
    plist = handles.zero_list(find(handles.pz == 1));
    plist = plist - j*imag(plist).*(1-handles.pairs(find(handles.pz == 1)));
    plist = [plist; conj(handles.zero_list(find(handles.pz == 1 & handles.pairs)))];
    
    [handles.B,handles.A] = zp2tf(zlist,plist,str2double(get(handles.Gain,'string')));
      
    if any(abs(plist) >= 1)
        set(handles.warning_text,'Visible','on');
    else
        set(handles.warning_text,'Visible','off');
    end
    
    [fr,w] = freqz(handles.B,handles.A,1024);
    %    w = w/2/pi*handles.fs;
    fr = abs(fr);
    
    if get(handles.LinearButton,'value')
        set(handles.fr,'YData',abs(fr));
    else
        warning off
        fr = 20*log10(fr);
        warning on
        %        m = mean(fr);
        %        fr(fr < m-50) = m-50;
        %        fr(fr > m+50) = m+50;
        set(handles.fr,'YData',fr);
    end        
    
    if handles.filter_eval == 2
        
        if length(handles.match_vals) == 1
            H = freqz(handles.B,handles.A,[handles.match_frq 1]);
            H = H(1);
        else
            H = freqz(handles.B,handles.A,(handles.match_frq*(1:length(handles.match_vals))));
        end
        lin_err = sqrt(mean((handles.match_vals - abs(H)).^2));
        warning off
        tmp = 20*log10(handles.match_vals) - 20*log10(abs(H));
        tmp(isinf(tmp) | isnan(tmp)) = [];
        db_err = sqrt(mean(tmp.^2));
        warning on
        
        set(handles.LinError,'string',num2str(lin_err));
        set(handles.dBError,'string',num2str(db_err));
        
    elseif handles.filter_eval == 1
        
        handles.passband = find(w < handles.cutoff(1));
        handles.stopband = find(w > handles.cutoff(2));
        
        handles.pass_min = min(fr(handles.passband));
        handles.pass_max = max(fr(handles.passband));
        handles.stop_max = max(fr(handles.stopband));
        
        axes(handles.MagResponse);
        ax = axis;
        set(handles.pmax_line,'YData',[handles.pass_max handles.pass_max ax(4)]);
        set(handles.pmin_line,'YData',[handles.pass_min handles.pass_min ax(3)]);
        set(handles.smin_line,'YData',[handles.stop_max handles.stop_max ax(4)]);
        
        if get(handles.LinearButton,'value')
            set(handles.LinError,'string',num2str(abs(20*log10(handles.pass_max/handles.pass_min))));
            set(handles.dBError,'string',num2str(20*log10(handles.pass_max/handles.stop_max)));
        else
            set(handles.LinError,'string',num2str(abs(handles.pass_max - handles.pass_min)));
            set(handles.dBError,'string',num2str(handles.pass_max - handles.stop_max));            
        end
    end    
    
    if get(handles.threeDon,'Value')
        
        ce = exp(-j*2*pi*(0:.001:1));
        cevals = abs(polyval(handles.B,ce)./polyval(handles.A,ce));
        
        vals = abs(polyval(handles.B(end:-1:1),1./handles.z)./polyval(handles.A(end:-1:1),1./handles.z));
        
        if get(handles.LinearButton,'value')
            vals(vals > max(cevals)+handles.upplus) = max(cevals)+handles.upplus;
            zmax = max([max(max(vals)),cevals+0.01]);
            zmin = 0;
        else
            warning off
            cevals = 20*log10(cevals);
            warning on
            m = mean(cevals(~isinf(cevals) & ~isnan(cevals)));
            cevals(cevals > 50+m) = m+50;
            cevals(cevals < m-50) = m-50;
            vals = 20*log10(abs(vals));
            vals(vals > max(cevals)+handles.log_upplus) = max(cevals)+handles.log_upplus;
            vals(vals < min(cevals)-handles.log_upplus) = min(cevals)-handles.log_upplus;
            zmax = max([max(max(vals)),cevals+1]);
            zmin = min([min(min(vals)),cevals-1]);
        end
        
        set(handles.ce1,'ZData',cevals + 0.01);
        set(handles.ce2,'ZData',cevals - 0.01);
        set(handles.surf,'ZData',vals);
        set(handles.surf,'CData',vals);        
        
    end
    
    guidata(h,handles);
    
    
elseif handles.Mode == 6
    
    a=handles.surf_plot;
    H=get(gcf,'currentpoint'); %  [ X     Y ]
    th(1)=atan(5*(H(1)-0.5)) - atan(5*(handles.G.pos(1)-0.5));
    th(2)=atan(5*(H(2)-0.5)) - atan(5*(handles.G.pos(2)-0.5));
    th=th*180/pi;
    newview=mod(handles.G.view-th,360);
    set(a,'view',newview);
    
end


% --------------------------------------------------------------------
function varargout = LinearButton_Callback(h, eventdata, handles, varargin)

if get(h,'value') == 1
    set(handles.DecibelButton,'value',0);
else
    set(handles.DecibelButton,'value',1);
end

if handles.filter_eval == 2
    if length(handles.match_vals) < 50
        set(handles.match_plot(1),'YData',handles.match_vals);
        set(handles.match_plot(2),'YData',kron(handles.match_vals,[0 1 NaN]));
    else
        set(handles.match_plot,'YData',handles.match_vals);
    end
end

axes(handles.surf_plot);
zlabel('|H(z)|');
axes(handles.MagResponse);
ylabel('|H(\omega)|');

slow_update_frq_resp(h,[],handles);



% --------------------------------------------------------------------
function varargout = DecibelButton_Callback(h, eventdata, handles, varargin)

if get(h,'value') == 1
    set(handles.LinearButton,'value',0);
else
    set(handles.LinearButton,'value',1);
end

if handles.filter_eval == 2
    if length(handles.match_vals) < 50
        mv = 20*log10(handles.match_vals);
        m = mean(mv);
        mv(mv < m-50) = m-50;
        mv(mv > m+50) = m+50;
        set(handles.match_plot,'YData',mv);
        set(handles.match_plot(1),'YData',mv);
        set(handles.match_plot(2),'YData',kron(mv,[0 1 NaN]));
    else
        mv = 20*log10(handles.match_vals);
        m = mean(mv);
        mv(mv < m-50) = m-50;
        mv(mv > m+50) = m+50;
        set(handles.match_plot,'YData',mv);
    end
end

axes(handles.surf_plot);
zlabel('|H(z)| (dB)');
axes(handles.MagResponse);
ylabel('|H(\omega)| (dB)');

slow_update_frq_resp(h,[],handles);



% --------------------------------------------------------------------
function varargout = FIRCoefs_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = LinError_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = dBError_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)

disp(' ');
disp('Filter coefficients exported to the variables ''A_pz'' and ''B_pz'':');

B_pz = handles.B;
while B_pz(end) == 0
    B_pz(end) = [];
end

A_pz = handles.A;
while A_pz(end) == 0
    A_pz(end) = [];
end

cstring = ['B_pz = [ ' num2str(B_pz,' %7g, ')];
cstring = deblank(cstring);
cstring(end:end+1) = ' ]';

disp(cstring);

cstring = ['A_pz = [ ' num2str(A_pz,'%7g, ')];
cstring = deblank(cstring);
cstring(end:end+1) = ' ]';

disp(cstring);

assignin('base','B_pz',B_pz);
assignin('base','A_pz',A_pz);

disp(' ');



% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)

disp(' ');
disp('Final filter coefficients (not exported to workspace):');


B_pz = handles.B;
while B_pz(end) == 0
    B_pz(end) = [];
end

A_pz = handles.A;
while A_pz(end) == 0
    A_pz(end) = [];
end

cstring = ['B_pz = [ ' num2str(B_pz,' %7g, ')];
cstring = deblank(cstring);
cstring(end:end+1) = ' ]';

disp(cstring);

cstring = ['A_pz = [ ' num2str(A_pz,'%7g, ')];
cstring = deblank(cstring);
cstring(end:end+1) = ' ]';

disp(cstring);

delete(gcf);


% --------------------------------------------------------------------
function varargout = PlaySound_Callback(h, eventdata, handles, varargin)

t = 0:(1/handles.fs):(2-1/handles.fs);
snd = 0;
for f=handles.match_frq:handles.match_frq:handles.fs/2
    snd = snd + cos(2*pi*t*f);
end

env = [linspace(0,1,handles.fs/8) ones(1,handles.fs*7/4) linspace(1,0,handles.fs/8)];
snd = snd.*env;

snd = filter(handles.B,handles.A,snd);

soundsc(snd,handles.fs);


% --------------------------------------------------------------------
function varargout = printgui_Callback(h, eventdata, handles, varargin)

printdlg -crossplatform


% --------------------------------------------------------------------
function varargout = ToClipboard_Callback(h, eventdata, handles, varargin)

print -dbitmap


% --------------------------------------------------------------------
function varargout = surf_plot_ButtonDownFcn(h, eventdata, handles, varargin)

handles.G.pos=get(gcf,'currentpoint'); % [ X    Y ]
handles.G.view=get(handles.surf_plot,'View');

handles.Mode = 6;
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = PZ_Place_3D_Figure_KeyPressFcn(h, eventdata, handles, varargin)

% right = 28
% left = 29
% up = 30
% down = 31

% move .005 per keypress

c = double(get(gcf,'CurrentCharacter'));

if handles.selected_zero > 0
    if c == 28 | c == 29 | c == 30 | c == 31
        z = handles.zero_list(handles.selected_zero);
        handles.Mode = 8;
        
        if c == 28
            handles.key_pt = [real(z)-.005 imag(z)];
        elseif c == 29
            handles.key_pt = [real(z)+.005 imag(z)];
        elseif c == 30 & handles.pairs(handles.selected_zero)
            handles.key_pt = [real(z) imag(z)+.005];
        elseif c == 31 & handles.pairs(handles.selected_zero)
            handles.key_pt = [real(z) imag(z)-.005];
        else
            handles.key_pt = [real(z) imag(z)];
        end
        
        figure1_WindowButtonMotionFcn(h,eventdata,handles);
    end
end




% --------------------------------------------------------------------
function varargout = threeDon_Callback(h, eventdata, handles, varargin)

if get(h,'Value')
    set(handles.surf,'Visible','On');
    set(handles.ce1,'Visible','On');
    set(handles.ce2,'Visible','On');
else
    set(handles.surf,'Visible','Off');
    set(handles.ce1,'Visible','Off');
    set(handles.ce2,'Visible','Off');
end

slow_update_frq_resp(h,[],handles);



% --------------------------------------------------------------------
function varargout = LpBp_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = HpBs_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = pair_cb_Callback(h, eventdata, handles, varargin)

sel = handles.selected_zero;
add_pz = 0;
if sel > 0
    if get(handles.pair_cb,'Value')
        handles.pairs(sel) = 1;
        set(handles.zero_handle_list(sel),'XData',real(handles.zero_list(sel))*[1 1],'YData',imag(handles.zero_list(sel))*[1 -1]);
        set(handles.ZeroMag,'string',num2str(abs(handles.zero_list(sel))));
        set(handles.ZeroAngle,'string',num2str(angle(handles.zero_list(sel))));
        add_pz = -2*handles.pz(sel)+1;
    else
        handles.pairs(handles.selected_zero) = 0;
        set(handles.zero_handle_list(sel),'XData',real(handles.zero_list(sel)),'YData',0);
        set(handles.ZeroMag,'string',num2str(abs(real(handles.zero_list(sel)))));
        set(handles.ZeroAngle,'string',num2str(angle(real(handles.zero_list(sel)))));
        add_pz = 2*handles.pz(sel)-1;
    end
else
    set(handles.pair_cb,'value',0);
end

if add_pz ~= 0
    handles.zero_list = [handles.zero_list; 0];
    handles.pairs = [handles.pairs; 0];
    handles.pz = [handles.pz; (add_pz + 1)/2];

    axes(handles.ZeroPlot);
    hdl2 = line([0],[0]);
    set(hdl2,'LineStyle','none','LineWidth',1,'MarkerSize',8,'Color',[0 1 0]);
    set(hdl2,'ButtonDownFcn',sprintf('pole_zero_place3d(''zero_select'',gcbo,[],guidata(gcbo));'));
    if handles.pz(end) == 0
        set(hdl2,'Marker','o');
    else
        set(hdl2,'Marker','x');
    end

    handles.zero_handle_list = [handles.zero_handle_list hdl2];
end

slow_update_frq_resp(h,[],handles);


% --- Executes on button press in AddPole.
function AddPole_Callback(hObject, eventdata, handles)
% hObject    handle to AddPole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AddPole
