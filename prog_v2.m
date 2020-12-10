clear all, clc
f = figure;
global gr;
set(f,'MenuBar','none')
set(f,'resize', "off")
set(f,'name', "Расчет потенциальной энергии взаимодействия атомов через волновые функции Шрёдингера")
set(f,'position',[300 200 1024 220]);

filem = uimenu ("label", "&File", "accelerator", "f");
editm = uimenu ("label", "&Actions", "accelerator", "e");
uimenu (filem, "label", "Export...", "accelerator", "q", "callback", "export_result()");
uimenu (filem, "label", "Exit", "accelerator", "q", "callback", "close (gcf)");
uimenu (editm, "label", "Start perebor...", "accelerator", "p", "callback", "start_perebor()");
#uimenu (editm, "label", "Start perebor2...", "accelerator", "p", "callback", "start_perebor2()");
uimenu (editm, "label", "Plot multiplot graph U(r)", "callback", "plotU_multiplot(1)");
uimenu (editm, "label", "Plot one graph U(r)", "callback", "plotU_multiplot(2)");
uimenu (editm, "label", "Plot graph F(r)", "callback", "plotFuse_multiplot()","separator", "on");
uimenu (editm, "label", "Plot animation graph U(r)", "callback", "export_animation_Ur()");
uimenu (editm, "label", "Plot points of graph U(r)", "callback", "export_animation_Ur_points_only()");

global MassOfAtom = containers.Map;
MassOfAtom("H")=1.00811;
MassOfAtom("He")=4.0026;
MassOfAtom("Li")=6.939;
MassOfAtom("Be")=9.0122;
MassOfAtom("B")=10.811;
MassOfAtom("C")=12.01115;
MassOfAtom("N")=14.0067;
MassOfAtom("O")=15.9994;
MassOfAtom("F")=18.9984;
MassOfAtom("Ne")=20.179;
global NumberOfAtom = containers.Map;
NumberOfAtom("H")=1;
NumberOfAtom("He")=2;
NumberOfAtom("Li")=3;
NumberOfAtom("Be")=4;
NumberOfAtom("B")=5;
NumberOfAtom("C")=6;
NumberOfAtom("N")=7;
NumberOfAtom("O")=8;
NumberOfAtom("F")=9;
NumberOfAtom("Ne")=10;

Electr_conf1 = "0-2-0-0";
Electr_conf2 = "0-2-0-0";
selected_atom1 = 1;
selected_atom2 = 1;
Xmin = 0; 
hx = 0.01;
Xmax = 30; 
%----------------------------
if version>="5.0.0" #сделано для совместимости c Octave 4.4.1
  res = isfile("value.txt");
else
  res = 1;%True
endif
if res
  load value.txt Electr_conf1 Electr_conf2 selected_atom1 selected_atom2 Xmin hx Xmax
endif  
%----------------------------

p1 = uipanel (f,"title", "For first atom:", "position", [0 0 .3 1]);
uicontrol(p1, "style", "text", "string", "Electr conf=", "position", [10 160 100 30]);
uicontrol(p1, "style", "text", "string", "Atom=", "position", [10 100 50 30]);
global Electr_conf1_g=uicontrol(p1, "style", "edit", "string", Electr_conf1, "position", [100 160 120 30]);
global Atom1_g=uicontrol(p1, "style", "popupmenu", "string", keys(NumberOfAtom), "position", [60 100 80 30],'tooltipstring', 'select first atom');

p2 = uipanel (f,"title", "For two atom:", "position", [0.3 0 .3 1]);
uicontrol(p2, "style", "text", "string", "Electr conf=", "position", [10 160 100 30]);
uicontrol(p2, "style", "text", "string", "Atom=", "position", [10 100 50 30]);
global Electr_conf2_g=uicontrol(p2, "style", "edit", "string", Electr_conf2, "position", [100 160 120 30]);
global Atom2_g=uicontrol(p2, "style", "popupmenu", "string", keys(NumberOfAtom), "position", [60 100 80 30],'tooltipstring', 'select second atom');

p3 = uipanel (f,"title", "Settings axis:", "position", [0.6 0 1 1]);
uicontrol(p3, "style", "text", "string", "Xmin=", "position", [10 160 40 30]);global Xmin_g=uicontrol(p3, "style", "edit", "string", num2str(Xmin), "position", [50 160 80 30]);
uicontrol(p3, "style", "text", "string", "Xmax=", "position", [140 160 40 30]);global Xmax_g=uicontrol(p3, "style", "edit", "string", num2str(Xmax), "position", [180 160 80 30]);
uicontrol(p3, "style", "text", "string", "hx=", "position", [260 160 40 30]);global hx_g=uicontrol(p3, "style", "edit", "string", num2str(hx), "position", [290 160 80 30]);
global CBenable_sigma_g = uicontrol(p3,'style','checkbox','position',[10 100 110 20],'string','sigma>0', 'tooltipstring', '1 -> sigma>0, 0 -> sigma=0');
global sigma_g_label = uicontrol(p3, "style", "text", "string", "sigma=", "position", [100 100 50 30]);
global sigma_g_value = uicontrol(p3,'style','edit','position',[150 100 80 25],'string','0', 'tooltipstring','value in Ang.');

global CBTypeUfunc_g = uicontrol(p3,'style','popupmenu','position',[10 130 150 20],'string',{'U4','U2','U1','U3: Molier^2 + voln func'}, 'tooltipstring', 'select U4 or U2');

global CBenable_hartree_fock_g = uicontrol(p3,'style','checkbox','position',[10 70 110 20],'string','Hartree-Fock', 'tooltipstring', 'for U4 only; Publication: ROOTHAAN-HARTREE-FOCK ATOMIC WAVEFUNCTION // ATOMIC DATA AND NUCLEAR DATA TABLES 14, 177-478 (1974)');
set(CBenable_hartree_fock_g,'Visible','off');


function selectionSigma(src,event)
  global sigma_g_label sigma_g_value Atom1_g Atom2_g MassOfAtom;      
  isChecked = get(src,'value');
  disp(['Selection: ' int2str(isChecked)]);
  if isChecked==false
    set(sigma_g_label,'Visible','off');
    set(sigma_g_value,'Visible','off');
    disp('hide')
  else
    set(sigma_g_label,'Visible','on');
    set(sigma_g_value,'Visible','on');
    [Electr_conf1, Electr_conf2, A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a0,a1,a2, selected_atom1, selected_atom2] = getMainValues()
    
    [sigma1, sigma2]=Get_Sigma(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a0,a1,a2, Atom1_g, Atom2_g,MassOfAtom, selected_atom1, selected_atom2)
    sigma=sqrt((sigma1^2+sigma2^2)/2);
    set(sigma_g_value,'string',num2str(sigma));
  endif
endfunction

function selectionhartree_fock(src,event)
  global CBTypeUfunc_g;      
  isChecked = get(src,'value');
  disp(['Selection: ' int2str(isChecked)]);
  if isChecked==false
    set(CBTypeUfunc_g,'Visible','on');
    disp('hide')
  else
    set(CBTypeUfunc_g,'Visible','off');
  endif
endfunction

set(CBenable_sigma_g,'callback', @selectionSigma);
selectionSigma(CBenable_sigma_g,0);

set(CBenable_hartree_fock_g,'callback', @selectionhartree_fock);
selectionhartree_fock(CBenable_hartree_fock_g,0);

global wbar = waitbar(0, 'Push plot for start','Visible', 'off');
c = get(wbar,'Children'); 
set(c,'Parent',p3, 'Position',[0.1 0.5 0.5 .05]);  % Set the position of the WAITBAR on your figure 

#-----load settings
#set(CBminimum_g,'value', type_minimum);
set(Atom1_g,'value', selected_atom1);
set(Atom2_g,'value', selected_atom2);

#-----------changes--------------
function [Electr_conf1, Electr_conf2, A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a0,a1,a2, selected_atom1, selected_atom2] = getMainValues()
  global Electr_conf1_g Electr_conf2_g Atom1_g Atom2_g NumberOfAtom
  printf("------------------------START---get main values----------------------\n");
  Electr_conf1 = get(Electr_conf1_g,'String');
  Electr_conf2 = get(Electr_conf2_g,'String');
  
  res1=sscanf(Electr_conf1,"%d-%d-%d-%d");
  A1 = res1(1);
  B1 = res1(2);
  C1 = res1(3);
  D1 = res1(4);
  res2=sscanf(Electr_conf2,"%d-%d-%d-%d");
  A2 = res2(1);
  B2 = res2(2);
  C2 = res2(3);
  D2 = res2(4);
  
  selected_atom1 = get(Atom1_g,'Value');
  selected_atom2 = get(Atom2_g,'Value');
  Z1 = NumberOfAtom(get(Atom1_g,'String'){selected_atom1});
  Z2 = NumberOfAtom(get(Atom2_g,'String'){selected_atom2});
    
  if sum(res1)!=Z1 || sum(res2)!=Z2  
    errordlg("Electronic configuration dont corrent! Please correct.","Error configuration structure of atoms");
    error("Error electronic configuration of atoms");
  endif
  
  a0 = 0.529;
  a1 = a0/Z1;
  a2 = a0/Z2;
  printf("------------------------END---get main values----------------------\n");
endfunction


function y = getSigmaddU(k,M)
  hc = 2000; #eV * Ang
  mc2 = M*1E+9; # eV
  y = sqrt(hc*sqrt(1.0/(k*mc2)));
endfunction 

function [sigma1, sigma2]=Get_Sigma(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a0,a1,a2, Atom1_g, Atom2_g,MassOfAtom, selected_atom1, selected_atom2)
  printf("------------------------START---Search sigma----------------------\n");
  [r_min2,y_min] = fminsearch(@(r) U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,a2,0),4)%второй минимум 
  
  k = ddU_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r_min2,a1,a2,0);
  printf("d^2U(r_min2)/dr^2=%f\n",k);
  M1 = MassOfAtom(get(Atom1_g,'String'){selected_atom1});
  M2 = MassOfAtom(get(Atom2_g,'String'){selected_atom2});
  sigma1= getSigmaddU(k,M1);
  sigma2= getSigmaddU(k,M2);
  printf("------------------------END---Search sigma----------------------\n");
endfunction



function y = U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,a2,sigma)
  global CBTypeUfunc_g CBenable_hartree_fock_g;
  TypeUfunc = get(CBTypeUfunc_g,'value');
  is_hartree_fock = get(CBenable_hartree_fock_g,'value');
  switch(TypeUfunc)
    case 1 #U4
      y = U4_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,a2,sigma,is_hartree_fock);    
    case 2 #U2
      y = U2_voln_func_sigma0(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1);  
    case 3#U1  
      y = U1_voln_func_sigma0(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1);
    case 4#U3 Molier + voln func
      y = U3_voln_func_Molier(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,a2);  
  end    
endfunction

function [x_arr,y_arr] = calcUOnly(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,sigma,a0, a1, a2, Xmax, Xmin, hx)
  global wbar
  x_arr=Xmin:hx:Xmax;
  step = 1/length(x_arr);
  for i=1:length(x_arr)
    y_arr(i)=U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x_arr(i)*a0,a1,a2,sigma);
    waitbar(step*i,wbar,sprintf('Progress %.2d %%. Please wait...',step*i*100));
  endfor
endfunction

function res = search_all_x(f, x, y)
  k1 = 1;
  for k=1:length(y)-1
    if abs(y(k))>1E-10
      #printf("yk=%f    yk+1=%f\n",y(k),y(k+1));
      if y(k)*y(k+1)<=0
          xk0(k1)=x(k);
          k1=k1+1;
      endif
    endif  
  endfor
  for k=1:length(xk0)
    xk(k) = fsolve(f,[xk0(k)],optimset('Display','off','TolX',1e-3));
  end
  res = xk;
endfunction

function plotU_multiplot(type)#построения графика потенциальной энергии взаимодействия атомов
  global Xmin_g hx_g Xmax_g Ymin_g Ymax_g wbar CBenable_sigma_g sigma_g_value;
  Xmin = str2double(get(Xmin_g,'String'));  Xmax = str2double(get(Xmax_g,'String'));  hx = str2double(get(hx_g,'String'));
  
  [Electr_conf1, Electr_conf2, A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a0,a1,a2, selected_atom1, selected_atom2] = getMainValues()
  
  enable_sigma = get(CBenable_sigma_g,'Value')
  if(enable_sigma)
    sigma=str2double(get(sigma_g_value,'string'))
  else
    sigma=0;
    sigma1=0;
    sigma2=0;
  endif
  printf("------------------------START Calc U(r)----------------------\n");  
  [x_arr,y_arr] = calcUOnly(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,sigma,a0, a1, a2, Xmax, Xmin, hx)
  waitbar(1,wbar,"Completed!");
  printf("------------------------END Calc U(r)----------------------\n");  
  printf("------------------------START Plot U(r)----------------------\n");  
  gr = figure('name','Graph of U(r)','Position', [10 100 1024 568])
  if type==1
    x_min_max = search_all_x(@(x) U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma), x_arr,y_arr)
    if length(x_min_max)==0
      error("Error: the function does not cross the OX axis");
    endif
    
    [xmin1,ymin1] = fminbnd(@(x) U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma),x_min_max(1),x_min_max(2)) %находим первый минимум
    [xmin2,ymin2] = fminbnd(@(x) U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma),x_min_max(3),Xmax) %находим второй минимум
    [xmax2,ymax2] = fminbnd(@(x) -U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma),x_min_max(2),x_min_max(3)) %находим первый максимум
    ymax2 = -ymax2;
    
    subplot(3,1,1);plot(x_arr,y_arr); xlim([0 x_min_max(2)]);ylim([ymin1 -ymin1]);xlabel('r/a_0');ylabel('U_4(r), eV');
    title(sprintf("a_0=%f Ang., x_{min2}=%f, U_{min2}= %f eV a_1=%f a_2=%f {\\sigma}=%f %s|%s",a0,xmin2,ymin2,a1,a2,sigma,Electr_conf1, Electr_conf2));
    subplot(3,1,2);plot(x_arr,y_arr); xlim([x_min_max(2) x_min_max(3)]);ylim([0 ymax2]);xlabel('r/a_0');ylabel('U_4(r), eV');
    subplot(3,1,3);plot(x_arr,y_arr); xlim([x_min_max(3)-1 xmin2+3]);ylim([ymin2 -ymin2]);xlabel('r/a_0');ylabel('U_4(r), eV');
  else#type==2
    [xmin1,ymin1] = fminbnd(@(x) U(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma),Xmin,Xmax) %находим первый минимум
    plot(x_arr,y_arr); xlim([Xmin Xmax]);xlabel('r/a_0');ylabel('U(r), eV');
    title(sprintf("a_0=%f Ang., x_{min1}=%f, U_{min1}= %f eV a_1=%f a_2=%f {\\sigma}=%f %s|%s",a0,xmin1,ymin1,a1,a2,sigma,Electr_conf1, Electr_conf2));
  endif  
  printf("------------------------END Plot U(r)----------------------\n");  
  save value.txt Electr_conf1 Electr_conf2 selected_atom1 selected_atom2 Xmin hx Xmax
endfunction

function plotFuse_multiplot()#построения графика силы F(r)=-dU/dr
  global Xmin_g hx_g Xmax_g Ymin_g Ymax_g wbar CBenable_sigma_g sigma_g_value;
  Xmin = str2double(get(Xmin_g,'String'));  Xmax = str2double(get(Xmax_g,'String'));  hx = str2double(get(hx_g,'String'));
  
  [Electr_conf1, Electr_conf2, A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a0,a1,a2, selected_atom1, selected_atom2] = getMainValues();
  
  enable_sigma = get(CBenable_sigma_g,'Value')
  if(enable_sigma)
    sigma=str2double(get(sigma_g_value,'string'));
  else
    sigma=0;
    sigma1=0;
    sigma2=0;
  endif
  printf("------------------------START Calc F(r)----------------------\n");  
  x_arr=Xmin:hx:Xmax;
  step = 1/length(x_arr);
  for i=1:length(x_arr)
    y_arr(i)=F_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x_arr(i)*a0,a1,a2,sigma);
    waitbar(step*i,wbar,sprintf('Progress %.2d %%. Please wait...',step*i*100));
  endfor
  waitbar(1,wbar,"Completed!");
  printf("------------------------END Calc F(r)----------------------\n");  
  printf("------------------------START Plot F(r)----------------------\n");  
  gr = figure('name','Graph of F(r)')
  x_min_max = search_all_x(@(x) F_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma), x_arr,y_arr)
  if length(x_min_max)==0
    error("Error: the function does not cross the OX axis");
  endif
  
  [xmin1,ymin1] = fminbnd(@(x) F_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma),x_min_max(1),x_min_max(2)) %находим первый минимум
  [xmin2,ymin2] = fminbnd(@(x) F_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma),x_min_max(3),Xmax) %находим второй минимум
  [xmax2,ymax2] = fminbnd(@(x) -F_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,x*a0,a1,a2,sigma),x_min_max(2),x_min_max(3)) %находим первый максимум
  ymax2 = -ymax2;
  
  subplot(3,1,1);plot(x_arr,y_arr); xlim([0 x_min_max(2)]);ylim([ymin1 -ymin1]);xlabel('r/a_0');ylabel('F(r), eV');
  title(sprintf("a_0=%f Ang., x_{min2}=%f, F_{min2}= %f eV a_1=%f a_2=%f sigma=%f",a0,xmin2,ymin2/a0,a1,a2,sigma));
  subplot(3,1,2);plot(x_arr,y_arr); xlim([x_min_max(2) x_min_max(3)]);ylim([-ymax2 ymax2]);xlabel('r/a_0');ylabel('F(r), eV');
  subplot(3,1,3);plot(x_arr,y_arr); xlim([x_min_max(3)-1 xmin2+3]);ylim([ymin2 -ymin2]);xlabel('r/a_0');ylabel('F(r), eV');
  printf("------------------------END Plot F(r)----------------------\n");  
  save value.txt Electr_conf1 Electr_conf2 selected_atom1 selected_atom2 Xmin hx Xmax sigma1 sigma2
endfunction


function export_result()
  try
    pkg load io
  catch
    msgbox(sprintf("Stop, because: %s",lasterr));
    error('IO package not loaded. Stop')
  end
  global Xmin_g hx_g Xmax_g Ymin_g Ymax_g CBenable_sigma_g wbar sigma_g_value;
  [file,path, filter] = uiputfile({'*.txt';'*.csv';'*.xlsx'},'File for save','data')
if isequal(file,0) || isequal(path,0)
   disp('User clicked Cancel.')
else
    disp(['User selected ',fullfile(path,file),' and then clicked Save.'])
    Xmin = str2double(get(Xmin_g,'String'));  Xmax = str2double(get(Xmax_g,'String'));  hx = str2double(get(hx_g,'String'));
    
    [Electr_conf1, Electr_conf2, A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a0,a1,a2, selected_atom1, selected_atom2] = getMainValues()
    
   enable_sigma = get(CBenable_sigma_g,'Value')
  if(enable_sigma)
    sigma=str2double(get(sigma_g_value,'string'));
  else
    sigma=0;
  endif
    
   [x_arr,y_arr] = calcUOnly(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,sigma,a0, a1, a2, Xmax, Xmin, hx);
   waitbar(1,wbar,"Completed!");
   if filter==2 %csv
      csvwrite(fullfile(path,file), [x_arr(:), y_arr(:)])
   elseif filter==3 %xlsx
      xlswrite(fullfile(path,file), [x_arr(:), y_arr(:)])
   else %txt
      dlmwrite(fullfile(path,file), [x_arr(:), y_arr(:)], "delimiter", " ", "newline", "\n");
   endif
   waitbar(1,wbar,"Save finished!");
end
endfunction

function start_perebor()#работает только для первых функций
  global Electr_conf1_g Electr_conf2_g Atom1_g Atom2_g CBenable_sigma_g NumberOfAtom wbar MassOfAtom;
  selected_atom1 = get(Atom1_g,'Value')
  selected_atom2 = get(Atom2_g,'Value')
  Z1 = NumberOfAtom(get(Atom1_g,'String'){selected_atom1})
  Z2 = NumberOfAtom(get(Atom2_g,'String'){selected_atom2})
  
  a0 = 0.529;
  aa1 = a0/Z1;
  aa2 = a0/Z2;
  i=1;
  step = 1/(Z1^2*Z2^2)^2;
  tic
  for a1=0:2
    for b1=0:2
      for c1=0:2
        for d1=0:4
          for a2=0:2
            for b2=0:2
              for c2=0:2
                for d2=0:4
                  if (a1+b1+c1+d1)==Z1 && (a2+b2+c2+d2)==Z2
                    [k_max,U_max] = fminsearch(@(k) -Uk_voln_func(k,a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,aa1,aa2,0),0.5);%второй минимум
                    str = sprintf("A1=%d B1=%d C1=%d D1=%d, A2=%d B2=%d C2=%d D2=%d, Z1=%d, Z2=%d, k_{max}=%f U_{kmax}=%f",a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,k_max,U_max)
                    result_Str{i} = str;
                    Umax_res(i) = U_max;
                    kmax_res(i) = k_max;
                    i=i+1;
                    if step<1
                      waitbar(step*i,wbar,sprintf('Found %.2d values. Please wait...',i));
                    endif  
                  endif
                endfor  
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor  
  endfor
  toc
  [M1,I1] = min(Umax_res)
  [M2,I2] = max(Umax_res)
  printf("1:kmax=%f U=%f\n",kmax_res(I1),M1);disp(result_Str{I1});
  printf("2:kmax=%f U=%f\n",kmax_res(I2),M2);disp(result_Str{I2});
  res=sscanf(result_Str{I1},"A1=%d B1=%d C1=%d D1=%d, A2=%d B2=%d C2=%d D2=%d, Z1=%d, Z2=%d, k_{max}=%f U_{kmax}=%f")
  set(Electr_conf1_g,'string',sprintf("%d-%d-%d-%d",res(1),res(2),res(3),res(4)));
  set(Electr_conf2_g,'string',sprintf("%d-%d-%d-%d",res(5),res(6),res(7),res(8)));
  answer = msgbox(sprintf("%d-%d-%d-%d | %d-%d-%d-%d; k_{max}=%f U_{kmax}=%f sigma=0",res(1),res(2),res(3),res(4),res(5),res(6),res(7),res(8),res(11),abs(res(12))),"Result");
  waitbar(1,wbar,"Completed!");
endfunction

function export_animation_Ur()
try
    pkg load io
  catch
    msgbox(sprintf("Stop, because: %s",lasterr));
    error('IO package not loaded. Stop')
end
try
    pkg load video
  catch
    msgbox(sprintf("Stop, because: %s",lasterr));
    error('Video package not loaded. Stop')
end
  global Xmin_g hx_g Xmax_g Atom1_g Atom2_g NumberOfAtom MassOfAtom CBenable_sigma_g gr wbar;
  [file_out,path, filter] = uiputfile({'*.gif';},'File for save','graph')
if isequal(path,0) || isequal(path,0)
   disp('User clicked Cancel.')
else
   disp(['User selected ',fullfile(path,file_out),' and then clicked Save.'])
   Xmin = str2double(get(Xmin_g,'String'))  
   Xmax = str2double(get(Xmax_g,'String'))  
   hx = str2double(get(hx_g,'String'))
   
   selected_atom1 = get(Atom1_g,'Value')
   selected_atom2 = get(Atom2_g,'Value')
   Z1 = NumberOfAtom(get(Atom1_g,'String'){selected_atom1})
   Z2 = NumberOfAtom(get(Atom2_g,'String'){selected_atom2})
   
   a0 = 0.529;
   aa1 = a0/Z1
   aa2 = a0/Z2
   enable_sigma = get(CBenable_sigma_g,'Value')
  
   fig = figure('name','График','DefaultAxesFontSize',10)
   waitbar(0,wbar,"Please wait...");
   tmp_file_pdf = fullfile(tempdir,"data.pdf")
   try
      delete(tmp_file_pdf)
   catch
      msgbox(sprintf("Dont delete file, because: %s",lasterr));
   end
   % цикл анимации
    for a1=0:2
      for b1=0:2
        for c1=0:2
          for d1=0:4
            for a2=0:2
              for b2=0:2
                for c2=0:2
                  for d2=0:4
                    if (a1+b1+c1+d1)==Z1 && (a2+b2+c2+d2)==Z2
                      pause(0.01);
                      if(enable_sigma)
                        [sigma1, sigma2]=Get_Sigma(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,a0,aa1,aa2, Atom1_g, Atom2_g,MassOfAtom, selected_atom1, selected_atom2)
                        sigma=sqrt((sigma1^2+sigma2^2)/2);
                      else
                        sigma=0;
                        sigma1=0;
                        sigma2=0;
                      endif
                      [x_arr,y_arr] = calcUOnly(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,sigma,a0, aa1, aa2, Xmax, Xmin, hx);#pbaspect([20 20 20]);
                      x_min_max = search_all_x(@(x) U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma), x_arr,y_arr)
                      if length(x_min_max)==0
                        error("Error: the function does not cross the OX axis");
                      endif
                      
                      [xmin1,ymin1] = fminbnd(@(x) U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma),x_min_max(1),x_min_max(2)) %находим первый минимум
                      [xmin2,ymin2] = fminbnd(@(x) U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma),x_min_max(3),Xmax) %находим второй минимум
                      [xmax2,ymax2] = fminbnd(@(x) -U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma),x_min_max(2),x_min_max(3)) %находим первый максимум
                      ymax2 = -ymax2;
                      
                      subplot(3,1,1);plot(x_arr,y_arr); xlim([0 x_min_max(2)]);ylim([ymin1 -ymin1]);xlabel('r/a_0');ylabel('U_4(r), eV');legend(sprintf("%d-%d-%d-%d | %d-%d-%d-%d",a1,b1,c1,d1,a2,b2,c2,d2));
                      title(sprintf("a_0=%f Ang., x_{min2}=%f, U_{min2}= %f eV a_1=%f a_2=%f {\\sigma}=%f",a0,xmin2,ymin2,aa1,aa2,sigma));
                      subplot(3,1,2);plot(x_arr,y_arr); xlim([x_min_max(2) x_min_max(3)]);ylim([0 ymax2]);xlabel('r/a_0');ylabel('U_4(r), eV');
                      subplot(3,1,3);plot(x_arr,y_arr); xlim([x_min_max(3)-1 xmin2+3]);ylim([ymin2 -ymin2]);xlabel('r/a_0');ylabel('U_4(r), eV');
                      drawnow
                      try
                        print(tmp_file_pdf,"-append")
                      catch
                        msgbox(sprintf("Stop, because: %s",lasterr));
                        error('error....stop...')
                      end
                   endif   
                endfor  
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor  
  endfor
  im = imread(tmp_file_pdf, "Index", "all");
  imwrite(im, fullfile(path,file_out), "DelayTime", .5)
  waitbar(1,wbar,"Completed!");
endif
endfunction
#############################################################################################################################################################################


function export_animation_Ur_points_only()
try
    pkg load io
  catch
    msgbox(sprintf("Stop, because: %s",lasterr));
    error('IO package not loaded. Stop')
end
try
    pkg load video
  catch
    msgbox(sprintf("Stop, because: %s",lasterr));
    error('Video package not loaded. Stop')
end
  global Xmin_g hx_g Xmax_g Atom1_g Atom2_g NumberOfAtom MassOfAtom CBenable_sigma_g gr wbar CBTypeUfunc_g;
  [file_out,path, filter] = uiputfile({'*.pdf';},'File for save','graph')
if isequal(path,0) || isequal(path,0)
   disp('User clicked Cancel.')
else
   disp(['User selected ',fullfile(path,file_out),' and then clicked Save.'])
   typeUfunc = get(CBTypeUfunc_g,'value');
   
   Xmin = str2double(get(Xmin_g,'String'))  
   Xmax = str2double(get(Xmax_g,'String'))  
   hx = str2double(get(hx_g,'String'))
   
   selected_atom1 = get(Atom1_g,'Value')
   selected_atom2 = get(Atom2_g,'Value')
   Z1 = NumberOfAtom(get(Atom1_g,'String'){selected_atom1})
   Z2 = NumberOfAtom(get(Atom2_g,'String'){selected_atom2})
   
   a0 = 0.529;
   aa1 = a0/Z1
   aa2 = a0/Z2
   enable_sigma = get(CBenable_sigma_g,'Value')
   i=1;
   waitbar(0,wbar,"Please wait...");
   % цикл анимации
    for a1=0:2
      for b1=0:2
        for c1=0:2
          for d1=0:4
            for a2=0:2
              for b2=0:2
                for c2=0:2
                  for d2=0:4
                    if (a1+b1+c1+d1)==Z1 && (a2+b2+c2+d2)==Z2
                      if(enable_sigma)
                        [sigma1, sigma2]=Get_Sigma(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,a0,aa1,aa2, Atom1_g, Atom2_g,MassOfAtom, selected_atom1, selected_atom2)
                        sigma=sqrt((sigma1^2+sigma2^2)/2);
                      else
                        sigma=0;
                        sigma1=0;
                        sigma2=0;
                      endif
                      
                      if(typeUfunc==1)#for U4
                        [x_arr,y_arr] = calcUOnly(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,sigma,a0, aa1, aa2, Xmax, Xmin, hx);#pbaspect([20 20 20]);
                        x_min_max = search_all_x(@(x) U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma), x_arr,y_arr)
                        if length(x_min_max)==0
                          error("Error: the function does not cross the OX axis");
                        endif
                      
                        [xmin1,ymin1] = fminbnd(@(x) U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma),x_min_max(1),x_min_max(2)) %находим первый минимум
                        [xmin2,ymin2] = fminbnd(@(x) U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma),x_min_max(3),Xmax) %находим второй минимум
                        [xmax2,ymax2] = fminbnd(@(x) -U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma),x_min_max(2),x_min_max(3)) %находим первый максимум
                        ymax2 = -ymax2;
                        
                        xmin_arr2(i) = xmin2;
                        ymin_arr2(i) = ymin2;
                        
                        xmax_arr1(i) = xmax2;
                        ymax_arr1(i) = ymax2;
                      else#for U2
                        [xmin1,ymin1] = fminbnd(@(x) U(a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,x*a0,aa1,aa2,sigma),Xmin,Xmax) %находим первый минимум
                      endif
                      
                      xmin_arr1(i) = xmin1;
                      ymin_arr1(i) = ymin1;
                      
                      El_conf{i} = sprintf("%d-%d-%d-%d | %d-%d-%d-%d",a1,b1,c1,d1,a2,b2,c2,d2);
                      i=i+1;
                   endif   
                endfor  
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor  
  endfor
  %--------------------------------------------------------------------------
  N = i-1
  disp("-------------------------")
  %cmap = hsv(N);
  fig = figure('name','График','DefaultAxesFontSize',10);
  disp("plot #1");
  if(typeUfunc==1)#for U4
      subplot(3,1,1);
      hold on;
      for j=1:N  
        if(xmin_arr2(j)>5.8 && xmin_arr2(j)<=6)
          marker = 'x';
          printf("%s: (%f;%f)\n",El_conf{j},xmin_arr1(j),ymin_arr1(j));
        elseif(xmin_arr2(j)>5.1 && xmin_arr2(j)<=5.6)
          marker = '>';
          printf("%s: (%f;%f)\n",El_conf{j},xmin_arr1(j),ymin_arr1(j));
        else
          marker = 'square';
        endif
        if(xmin_arr1(j)>2.4 && xmin_arr1(j)<=2.85 && ymin_arr1(j)>-1.5)% && ymin_arr1(j)>=-1.5
          printf("found: %s (%f;%f)\n",El_conf{j},xmin_arr1(j),ymin_arr1(j));
        endif    
        plot(xmin_arr1(j),ymin_arr1(j),'Marker',marker,'linestyle','none','Color','blue');
      endfor
      hold off;
      xlabel('r/a_0');ylabel('U_4(r), eV');
      
      disp("plot #2");
      subplot(3,1,2);
      hold on;
      for j=1:N  
        if(xmin_arr2(j)>5.8 && xmin_arr2(j)<=6)
          marker = 'x';
          printf("%s: (%f;%f)\n",El_conf{j},xmax_arr1(j),ymax_arr1(j));
        elseif(xmin_arr2(j)>5.1 && xmin_arr2(j)<=5.6)
          marker = '>';
          printf("%s: (%f;%f)\n",El_conf{j},xmax_arr1(j),ymax_arr1(j));
        else
          marker = 'square';
        endif  
        plot(xmax_arr1(j),ymax_arr1(j),'Marker',marker,'linestyle','none','Color','green');%
      endfor
      hold off;
      xlabel('r/a_0');ylabel('U_4(r), eV');
      
      disp("plot #3");
      subplot(3,1,3);
      hold on;
      for j=1:N  
        if(xmin_arr2(j)>5.8 && xmin_arr2(j)<=6)
          printf("%s: (%f;%f)\n",El_conf{j},xmin_arr2(j),ymin_arr2(j));
          marker = 'x';
        elseif(xmin_arr2(j)>5.1 && xmin_arr2(j)<=5.6)
          marker = '>';
          printf("%s: (%f;%f)\n",El_conf{j},xmin_arr2(j),ymin_arr2(j));
        else
          marker = 'square';
        endif  
        plot(xmin_arr2(j),ymin_arr2(j),'Marker',marker,'linestyle','none','Color','red');
      endfor
      hold off;
      xlabel('r/a_0');ylabel('U_4(r), eV');
      dlmwrite(fullfile(path,"data_points_extremum.txt"), [xmin_arr1(:), ymin_arr1(:), xmax_arr1(:),ymax_arr1(:), xmin_arr2(:), ymin_arr2(:)], "delimiter", " ", "newline", "\n");
  else
    plot(xmin_arr1,ymin_arr1,'Marker','square','linestyle','none','Color','blue');
    xlabel('r/a_0');ylabel('U_2(r), eV');    
    dlmwrite(fullfile(path,"data_points_extremum.txt"), [xmin_arr1(:), ymin_arr1(:)], "delimiter", " ", "newline", "\n");
  endif 
  print(fullfile(path,file_out));
  waitbar(1,wbar,"Completed!");
endif
endfunction
##############################################################################################################################################################################
function y=Uk2(k,A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a1,a2,sigma)
  F1=@(k) 16*A1./(4+a1.^2.*k.^2).^2.+B1*(1.0-3*a1.^2.*k.^2+2*a1.^4.*k.^4)./(a1.^2.*k.^2+1.0).^4 .+C1*(1.0-5*a1.^2.*k.^2)./(a1.^2.*k.^2+1.0).^4.+D1./(1.0+a1.^2.*k.^2).^3;
  F2=@(k) 16*A2./(4+a2.^2.*k.^2).^2.+B2.*(1.0-3*a2.^2.*k.^2+2*a2.^4.*k.^4)./(a2.^2.*k.^2+1.0).^4 .+C2.*(1.0-5*a2.^2.*k.^2)./(a2.^2.*k.^2+1.0).^4.+D2./(1.0+a2.^2.*k.^2).^3;
  y = (1.0-F1(k)).^2.*(1-F2(k)).^2.*exp(-k.^2.*sigma.^2)./k;
endfunction

function res = solve_dihotomi(f, a, b, eps)
 iter=0;
 while (abs(b-a)>eps && iter<20000)
   c=(a + b)/2;
   if(f(a)*f(c)<=0) 
      b = c;
   else
      a = c;
   endif
   iter=iter+1;
 end
 res=c;
endfunction

function start_perebor2()
  global Electr_conf1_g Electr_conf2_g Atom1_g Atom2_g CBenable_sigma_g NumberOfAtom wbar MassOfAtom;
  selected_atom1 = get(Atom1_g,'Value')
  selected_atom2 = get(Atom2_g,'Value')
  Z1 = NumberOfAtom(get(Atom1_g,'String'){selected_atom1})
  Z2 = NumberOfAtom(get(Atom2_g,'String'){selected_atom2})
  
  a0 = 0.529;
  aa1 = a0/Z1;
  aa2 = a0/Z2;
  i=1;
  h=0.1;
  tic
  for a1=0:h:1
    for b1=0:h:1
      for c1=0:h:1
        for d1=0:h:1
          for a2=0:h:1
            for b2=0:h:1
              for c2=0:h:1
                for d2=0:h:1
                  if (a1+b1+c1+d1)==1 && (a2+b2+c2+d2)==1
                    #[k_max,U_max] = fminsearch(@(k) -Uk2(k,a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,aa1,aa2,0),0.5);%второй минимум
                    param(1) = a1;
                    param(2) = b1;
                    param(3) = c1;
                    param(4) = d1;
                    param(5) = a2;
                    param(6) = b2;
                    param(7) = c2;
                    param(8) = d2;
                    param(9) = Z1;
                    param(10) = Z2;
                    param(12) = aa1;
                    param(13) = aa2;
                    f = @(k) (Uk2(k+0.001,a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,aa1,aa2,0)-Uk2(k-0.001,a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,aa1,aa2,0))/(2.0*0.001);
                    k_max = solve_dihotomi(f, 0.1, 100, 1E-5);
                    U_max = -Uk2(k_max,a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,aa1,aa2,0);
                    str = sprintf("A1=%f B1=%f C1=%f D1=%f, A2=%f B2=%f C2=%f D2=%f, Z1=%d, Z2=%d, k_{max}=%f U_{kmax}=%f",a1,b1,c1,d1,a2,b2,c2,d2,Z1,Z2,k_max,U_max)
                    result_Str{i} = str;
                    Umax_res(i) = U_max;
                    kmax_res(i) = k_max;
                    i=i+1;
                    waitbar(0,wbar,sprintf('Found %.2d values. Please wait...',i));
                  endif
                endfor  
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor  
  endfor
  toc
  [M,I] = min(Umax_res)
  printf("kmax=%f U=%f\n",kmax_res(I),M)
  disp(result_Str{I})
  res=sscanf(result_Str{I},"A1=%f B1=%f C1=%f D1=%f, A2=%f B2=%f C2=%f D2=%f, Z1=%d, Z2=%d, k_{max}=%f U_{kmax}=%f")
  set(Electr_conf1_g,'string',sprintf("%.2f-%.2f-%.2f-%.2f",res(1),res(2),res(3),res(4)));
  set(Electr_conf2_g,'string',sprintf("%.2f-%.2f-%.2f-%.2f",res(5),res(6),res(7),res(8)));
  answer = msgbox(sprintf("%.2f-%.2f-%.2f-%.2f | %.2f-%.2f-%.2f-%.2f; k_{max}=%f U_{kmax}=%f",res(1),res(2),res(3),res(4),res(5),res(6),res(7),res(8),res(11),abs(res(12))),"Result");
  waitbar(1,wbar,"Completed!");
endfunction