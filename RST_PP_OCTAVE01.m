function RST_PP_OCTAVE01
  pkg load image;
  hwnd=figure;
  %clear hwnd;

  dx_large=1/3; dy_large=1/3;
  dx=dx_large/6; dy=dy_large/7;
  dx_edit=dx*4;
  
  %MAKE GUI CONTROLS
  xpos=0; ypos=0;%column 1
  hTs=makeField('Ts', xpos,ypos,dx,dx,dy);
  hDelay=makeField('Delay', xpos+dx*3,ypos,dx,dx,dy);
  ypos+=dy;
  hNum=makeField('Plant num.', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hDen=makeField('Plant den.', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hContPlant=makeCheckbox('Continuous plant', [xpos+dx ypos dx_edit dy]);
  ypos+=dy;
  hBp=makeField('B(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hAp=makeField('A(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hUpdate=makeButton('Update plant only', [xpos+dx ypos dx_edit dy], @update_callback);
  
  xpos=dx_large; ypos=0;%column 2
  hRegW0=makeField('Reg. w0', xpos,ypos,dx,dx,dy);
  hRegZeta=makeField('Reg. zeta', xpos+dx*3,ypos,dx,dx,dy);
  ypos+=dy;
  hTrackW0=makeField('Track. w0', xpos,ypos,dx,dx,dy);
  hTrackZeta=makeField('Track. zeta', xpos+dx*3,ypos,dx,dx,dy);
  ypos+=dy;
  hContDesign=makeCheckbox('Continuous design criteria', [xpos+dx ypos dx_edit dy]);
  ypos+=dy;
  hP=makeField('P(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hBm=makeField('Bm(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hAm=makeField('Am(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hSolve=makeButton('Solve', [xpos+dx ypos dx_edit dy], @solve_callback);
  
  xpos=dx_large*2; ypos=0;%column 3
  hRemoveZeros=makeCheckbox('Remove stable zeros', [xpos ypos dx*2.5 dy]);
  hContCurve=makeCheckbox('Continuous curve', [xpos+dx*2.5 ypos dx*2.5 dy]);
  ypos+=dy;
  hHr=makeField('Hr(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hHs=makeField('Hs(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hDistFreq=makeField('Dist. freq.', xpos,ypos,dx,dx,dy);
  ypos+=dy;
  hR=makeField('R(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hS=makeField('S(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  hT=makeField('T(z)', xpos,ypos,dx,dx_edit,dy);
  ypos+=dy;
  
  xpos=0; ypos=dy_large;
  hOutput=makeAxes([xpos+dx ypos+dy 1-dx*2 dy_large-dy*2]);
  ypos+=dy_large;
  hControlSignal=makeAxes([xpos+dx ypos+dy 1-dx*2 dy_large-dy*2]);
  
  %SET DEFAULT VALUES
  set(hTs, 'string','1');
  set(hContPlant, 'value',0);%select discrete plant
  
    %continuous plant
  set(hDelay, 'string','0.6');
  set(hNum, 'string','1');
  set(hDen, 'string','[10 1]');
  
    %discrete plant
  set(hBp, 'string','[0 0.1 0.2]');
  set(hAp, 'string','[1 -1.3 0.42]');
  
  set(hContDesign, 'value',1);%check 'continuous design criteria'
  
    %regulation P(z)
  set(hRegW0, 'string','0.4');
  set(hRegZeta, 'string','0.9');
  set(hP, 'string','[1 -1.3741 0.4867]');
  
    %tracking Bm(z)/Am(z)
  set(hTrackW0, 'string','0.5');
  set(hTrackZeta, 'string','0.9');
  set(hBm, 'string','[0.0927 0.0687]');
  set(hAm, 'string','[1 -1.2451 0.4066]');
  
  set(hHr, 'string','1');
  set(hHs, 'string','1');
  set(hDistFreq, 'string','0');%disturbance cosine frequency
  
  set(hRemoveZeros, 'value',0);%don't cancel stable zeros
  
  solve_callback();
  
  
  
  
  function handle=makeText(text, pos)
    handle=uicontrol(hwnd, 'style','text', 'string',text, 'fontsize',8,...
      'units','normalized', 'position',[pos(1) 1-pos(2)-pos(4) pos(3) pos(4)]);
  endfunction
  function handle=makeEdit(pos)
    handle=uicontrol(hwnd, 'style','edit',...
      'units','normalized', 'position',[pos(1) 1-pos(2)-pos(4) pos(3) pos(4)]);
  endfunction
  function hEdit=makeField(name, x, y, textWidth, editWidth, height)
    hText=makeText(name, [x y textWidth height]);
    hEdit=makeEdit([x+textWidth y editWidth height]);
  endfunction
  
  function handle=makeCheckbox(name, pos)
    handle=uicontrol(hwnd, 'style','checkbox', 'string',name, 'fontsize',8,...
      'units','normalized', 'position',[pos(1) 1-pos(2)-pos(4) pos(3) pos(4)]);
  endfunction
  function handle=makeButton(name, pos, cb)
    handle=uicontrol(hwnd, 'string',name, 'fontsize',8, 'callback',cb,...
      'units','normalized', 'position',[pos(1) 1-pos(2)-pos(4) pos(3) pos(4)]);
  endfunction
  function handle=makeAxes(pos)
    handle=axes(hwnd, 'units','normalized', 'position',[pos(1) 1-pos(2)-pos(4) pos(3) pos(4)]);
  endfunction
  
  
  
  
  function P=z_polynomial(w0, zeta, Ts)
    s=w0*(-zeta+1i*sqrt(1-zeta*zeta));
    z=exp(s*Ts);
    rez=real(z); imz=imag(z);
    P=[1, -2*rez, rez*rez+imz*imz];
  endfunction

  function lz=count_leading_zeros(v)
    size_v=size(v);
    if size_v(1)<size_v(2)%horizontal vector
        v=v';
    end
    size_v=size(v);
    lz=0;
    for k=1:size_v(1)
        if v(k)==0
            lz=lz+1;
        else
            break;
        end
    end
  endfunction

  function [Ts B A Hr Hs dist P Bm Am]=acquire_data()
    c=clock;
    iteration_time=sprintf('%d-%d-%d %d:%d:%f', c(1), c(2), c(3), c(4), c(5), c(6))

    Ts=str2double(get(hTs, 'string'));
    BW=0;%rad/s
    %GET PLANT PTF
    if get(hContPlant, 'value')%continuous plant
        num=str2num(get(hNum, 'string'));
        den=str2num(get(hDen, 'string'));
        delay=str2double(get(hDelay, 'string'));
        sysc=tf(num, den, 'iodelay', delay);
        sysd=c2d(sysc, Ts, 'zoh');
        [B A]=tfdata(sysd, 'v');
        
        %get bandwidth
        BW=max(abs(roots(den)));
        
        B=padarray(B', sysd.ioDelay, 0, 'pre')';
        %get discretized delay
        %r1=step(sysd, 0:Ts:100);
        %totaldelay=count_leading_zeros(r1);
        %Bdelay=count_leading_zeros(B);
        %B=padarray(B', totaldelay-Bdelay, 0, 'pre')';
        set(hBp, 'string',['[' num2str(B) ']']);
        set(hAp, 'string',['[' num2str(A) ']']);
    else                                 %discrete plant
        B=str2num(get(hBp, 'string'));
        A=str2num(get(hAp, 'string'));
    end
    %GET REGULATION P(z) & TRACKING Bm(z)/Am(z)
    if get(hContDesign, 'value')%continuous design criteria
      w0  =str2double(get(hRegW0, 'string'));
      zeta=str2double(get(hRegZeta, 'string'));
      P=z_polynomial(w0, zeta, Ts);
      set(hP, 'string',['[' num2str(P) ']']);
      BW=max([BW w0]);
      
      w0  =str2double(get(hTrackW0, 'string'));
      zeta=str2double(get(hTrackZeta, 'string'));
      Am=z_polynomial(w0, zeta, Ts);
      set(hAm, 'string',['[' num2str(Am) ']']);
      BW=max([BW w0]);
      
      Bm=sum(Am);%Bm(z)=Am(1)
      set(hBm, 'string',['[' num2str(Bm) ']']);
    else                            %discrete design criteria
      P =str2num(get(hP, 'string'));
      Am=str2num(get(hAm, 'string'));
      Bm=str2num(get(hBm, 'string'));
    end
    set(hwnd, 'NumberTitle','off', 'Name',sprintf('RST Pole Placement v2. Nyquist period=%fs.', pi/BW));%%

    Hr=str2num(get(hHr, 'string'));
    Hs=str2num(get(hHs, 'string'));
    distf=str2double(get(hDistFreq, 'string'));
    if distf==0
        dist=[1 -1];
    else
        dist=[1, -2*cos(2*pi*distf*Ts), 1];
    end
  endfunction

  function no_int=has_no_integrator(A)%finds if A has no integrator
    r=roots(A);
    sr=size(r);
    no_int=true;
    for k=1:sr(1)
        if r(k)==1
            no_int=false;
        end
    end
  endfunction

  function [B_stable B_unstable]=separate_B(B)
    lz=count_leading_zeros(B);
    gain=B(lz+1);%find polynomial gain
    r=roots(B);
    B_stable=1; B_unstable=1;
    sr=size(r);
    for k=1:sr(1)
        if abs(r(k))>=1 %unstable & critical zeros
            B_unstable=conv(B_unstable, [1 -r(k)]);
        else            %stable zeros
            B_stable=conv(B_stable, [1 -r(k)]);
        end
    end
    B_stable=B_stable*gain;%put the gain in the cancelled part
    B_unstable=padarray(B_unstable', lz, 0, 'pre')';%put the delays in non-cancelled part
  endfunction

  function [t X ref_size]=sim_response(sim_time, t_ref, t_dist, TF_ref, TF_dist, Ts)
    res_ref =step(TF_ref, 0:(sim_time-t_ref));
    res_dist=impulse(TF_dist, 0:(sim_time-t_dist));
    ref_size=size(res_ref);
    res_ref  =padarray(res_ref,  t_ref,  0, 'pre');
    res_dist =padarray(res_dist, t_dist, 0, 'pre');
    t=0:Ts:(sim_time*Ts);
    X=res_ref-0.25*res_dist;%disturbance is a negative quarter-step
  endfunction

  function simulate_RST(Ts, B, A, R, S, T, Bm, Am, dist)
    %CALCULATE TFs:
    track=filt(Bm, Am);
    AS=conv(A, S);
    BR=conv(B, R);
    size_AS=size(AS);
    size_BR=size(BR);
    if size_AS(2)>size_BR(2)%AS longer than BR
        BR=padarray(BR', size_AS(2)-size_BR(2), 0, 'post')';
    else
        if size_AS(2)<size_BR(2)%BR longer than AS
            AS=padarray(AS', size_BR(2)-size_AS(2), 0, 'post')';
        end
    end
    den=AS+BR;

    %'c(z)/r(z):'
    c_r=minreal(track*filt(conv(T, B), den));

    %'Output sensitivity fn c(z)/D(z):'
    c_D=minreal(filt(AS, den));

    %'U(z)/r(z):'
    U_r=minreal(track*filt(conv(T, A), den));

    %'U(z)/D(z):'
    U_D=minreal(filt(conv(A, -R), den));

    D=filt(1, dist);

    %RESPONSE PLOTS
    sim_time=100; t_ref=5; t_dist=50;

    continuous_curve=get(hContCurve, 'value');
    [t PlantOutput ref_size]=sim_response(sim_time, t_ref, t_dist, c_r, c_D*D, Ts);
    reference=ones(ref_size);
    reference=padarray(reference, t_ref, 0, 'pre');
    %stairs(hOutput, t,PlantOutput);
    if continuous_curve
        plot(hOutput, t,reference,'-', t,PlantOutput,'-');
    else
        plot(hOutput, t,reference,'.', t,PlantOutput,'o');
    end
    set(hOutput, 'XMinorGrid','on', 'YMinorGrid','on');
    title(hOutput, sprintf('Plant output, step disturbance at %g s', t_dist*Ts));

    [t ControlSignal ref_size]=sim_response(sim_time, t_ref, t_dist, U_r, U_D*D, Ts);
    %stairs(hControlSignal, t,ControlSignal);
    if continuous_curve
        plot(hControlSignal, t,ControlSignal,'-');
    else
        plot(hControlSignal, t,ControlSignal,'o');
    end
    set(hControlSignal, 'XMinorGrid','on', 'YMinorGrid','on');
    title(hControlSignal, sprintf('Control signal, step disturbance at %g s', t_dist*Ts));
  endfunction

  function solve_callback()
    %'solve'
    % {
    [Ts Bp Ap Hr Hs dist P Bm Am]=acquire_data();

    Bps=[]; Bpu=[];
    if get(hRemoveZeros, 'value')%cancel stable zeros
        [Bps Bpu]=separate_B(Bp);
    else                         %don't cancel
        Bps=1; Bpu=Bp;
    end

    %SOLUTION
    %Ap Hs S + Bp Hr R = P      (no zeros cancelled)
    S_has_integrator=has_no_integrator(Ap);
    Hs_full=[];
    if S_has_integrator
        Hs_full=conv(Hs, [1 -1]);
    else
        Hs_full=Hs;
    end
    A=conv(Ap, Hs_full)';
    B=conv(Bpu, Hr)';
    nA=size(A); nA=nA(1)-1;%vertical vector size, power of z = size-1
    nB=size(B); nB=nB(1)-1;
    A=padarray(A, nB-1, 0, 'post');
    B=padarray(B, nA-1, 0, 'post');
    M=[];
    for k=1:nB
        M=[M A];
        A=circshift(A, 1);
    end
    for k=(nB+1):(nB+nA)
        M=[M B];
        B=circshift(B, 1);
    end
    %M
    n=size(M);  n=n(1);%M is a square matrix
    nP=size(P); nP=nP(2);
    P_full=padarray(P', n-nP, 0, 'post');
    SR=M\P_full;

    %Extract R, S, T
    R=SR((nB+1):(nB+nA))';
    S=SR(1:nB)';
    if S_has_integrator
        S=conv(S, [1, -1]);
    end
    S=conv(S, Bps);
    T=P/sum(Bpu);%(no zeros cancelled)
    set(hR, 'string',['[' num2str(R) ']']);%full R S T values for Update plant & Simulink buttons
    set(hS, 'string',['[' num2str(S) ']']);
    set(hT, 'string',['[' num2str(T) ']']);

    R=conv(R, Hr);
    S=conv(S, Hs);
    simulate_RST(Ts, Bp, Ap, R, S, T, Bm, Am, dist);
    %}
  endfunction
  function update_callback()
    [Ts Bp Ap Hr Hs dist P Bm Am]=acquire_data(handles);

    R=str2num(get(hR, 'string'));%full R S T values
    S=str2num(get(hS, 'string'));
    T=str2num(get(hT, 'string'));

    R=conv(R, Hr);
    S=conv(S, Hs);
    simulate_RST(Ts, Bp, Ap, R, S, T, Bm, Am, dist, handles);
  endfunction
endfunction
