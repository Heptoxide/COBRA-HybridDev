function dudx=derivative_dwt(u,wt name,wt level,dx,trt flag)
 %
 % Differentiation (Derivative) of Sampled Data Based on Discrete Wavelet ...
Transform

 %
 % dudx=derivative dwt(u,wt name,wt level,dx)
 %
 % u: uniformly?sampled data
 % wt name: name of the wavelet function (haar or spl)
 % wt level: level of wavelet decomposition
 % dx: sampling interval (default=1)
 % trt flag: flag of translation?rotation transformation for boundary ...
effect (default=1)
 % dudx: differentiations of data (u)
 %
 %See also
 % derivative cwt
 %
 % Reference:
 % J. W. Luo, J. Bai, and J. H. Shao,
 % "Application of the wavelet transforms on axial strain calculation in ...
ultrasound elastography,"
 % Prog. Nat. Sci., vol. 16, no. 9, pp. 942?947, 2006.

 if nargin<5
 trt flag=1;
 end
 if nargin<4
 dx=1;
 end

 if trt flag
 x=(1:length(u))*dx;
 a=(u(end)?u(1))/(x(end)?x(1));
 b=u(1)?a*x(1);
 u=u?a*x?b;
 else
 a=0;
 end


 wt name=lower(wt name);

 if strcmp(wt name,'haar')
 h0=[sqrt(2)/2 sqrt(2)/2]; %the decomposition low?pass filter
 h1=[?sqrt(2)/2 sqrt(2)/2]; %the decomposition high?pass filter
 elseif strcmp(wt name,'spl')
 h0=[0.125 0.375 0.375 0.125]*sqrt(2);
 h1=[?2 2]*sqrt(2);
 else
 error('wavelet name error');
 end

 y0=u;

 % Algorithme a Trous
 for n=1:wt level
 h0 atrous=[h0' zeros(length(h0),2ˆ(n?1)?1)]';
 h0 atrous=h0 atrous(1:(length(h0)?1)*(2ˆ(n?1)?1)+length(h0));

 h1 atrous=[h1' zeros(length(h1),2ˆ(n?1)?1)]';
 h1 atrous=h1 atrous(1:(length(h1)?1)*(2ˆ(n?1)?1)+length(h1));

 y1=conv(y0,h1 atrous);
 y0=conv(y0,h0 atrous);
 end

 index=round(length(y1)/2?length(u)/2)+[1:length(u)];
 dudx=y1(index);

 wt scale=2ˆwt level;

 if strcmp(wt name,'haar')
 dudx=?dudx/wt scaleˆ(3/2)*4;
 elseif strcmp(wt name,'spl')

 dudx=?dudx/wt scaleˆ(3/2);
 else
 error('wavelet name error');
 end

 dudx=dudx/dx+a;