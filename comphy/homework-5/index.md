---
mode: passage
title: 计算物理第五次作业
date: 2023/1/10
author: gky
# header: true
---
[作业](./HW5.pdf)
[解答](./%E8%A7%A3%E7%AD%94.pdf)
# 计算物理第五次作业
## 1.数值求解KdV方程
对KdV方程$u_t+uu_x+\del^2u_{xxx} = 0$作离散化处理：
$$
\frac{u_i^{j+1} -u_i^{j-1}}{2\D t} +\frac 1 3 (u_{i-1}^j+u_i^j+u_{i+1}^j)\frac{u_{i+1}^j-u_{i-1}^j}{2\D x}+\frac {\del^2(u_{i+2}^j-2u_{i+1}^j+2u_i^j-u_{i-1}^j)}{2\D x^3}=0
$$
其中时间和空间项均采用中心差分处理，并把$u$处理为相邻三项的平均值。
进行数值模拟，得到的结果如[视频](T1_delta=0.042,dt=0.001,nx=100.mp4)所示:
KdV方程在给定初始状态下演化时，会先产生一个较大的落差，之后在其两侧产生波动，并逐渐变成几个波峰同时存在的状态：
\img["T1_0.28s.png"]
\img["T1_0.613s.png"]
\img["T1_1.256s.png"]
之后，几个波包各自移动，相互不发生干扰：下图中显示了$t=2.043s$的情形，可以看到，最高的波包和最低的相互穿过对方，而未产生其他影响。
\img["T1_2.043s.png"]
在演化时间范围内，图像有一定抖动，但总体稳定性良好
## 2.数值求解扩散方程

对T(x,y,t)在空间坐标上做傅里叶变换,将方程转变为$\tilde{T}(\v k,t)$关于$t$的常微分方程,求解后作逆傅里叶变换,得到解的动画如视频[T2.mp4](./T2.mp4)所示.视频中取初态为以中心$(0.5,0.5)$为原点的高斯分布.
由于在上面的热力图中很难看出当t很大时的行为,故通过让colormap范围随着区域数值的最大值进行变动,得到新的动画[T2_.mp4](T2_.mp4),可以很好地反映长时间温度分布的状况
可以看到,在演化过程中,函数逐渐向左右两侧"扩散",最终形成了横向条状的温度分布.这是因为上下两端是绝热的(${\partial T\over\partial y}|_{y=0,1} =0$),而左右两端与低温热源接触($T|_{x=0,1}=0$),最终导致热量向左右两侧流失.

## 3.波动方程的求解
#### (1)
由波动方程及其边界条件$u|_{\text{boundary}}=0,0\le x,y\le 1$可利用本征函数展开得到通解:
$$
u(x,y,t) = \sum_{m,n=1}^{\infty } A_{mn}(t)\sin m\pi x\sin n\pi y\quad x,y\in [0,1],t\ge 0
$$
代入初值条件:
$$
A_{mn}(0)=\begin{cases}
1, & m=1,n=2\\
0, & \text{otherwise}
\end{cases}\\
A_{mn}'(0) = 0
$$
各个参数满足常微分方程:
$$
A''_{mn}(t)= \lambda (-m^2-n^2)\pi^2 A_{mn}(t)
$$
解出$A_{12}(t)=\cos\sqrt{5\lambda}\pi t;A_{mn}(t)=0,\text{otherwise}$
$\lambda = 1$得到
$$u(x,y,t) = \cos\sqrt{5 }\pi t\sin \pi x\sin 2\pi y
$$
#### (2)
对方程作离散化处理:定义$u_{ij}^k \equiv u(i\D x,j\D y,k\D t),0\le i< N_x,0\le i< N_y,0\le t<N_t$
$$
\frac{u_{ij}^{k+1}-2u_{ij}^k+u_{ij}^{k-1}}{(\D t)^2} = \lambda(\frac{u_{i+1,j}^k-2u_{ij}^k+u_{i-1,j}^k}{(\D x)^2}+\frac{u_{i,j+1}^k-2u_{ij}^k+u_{i,j-1}^k}{(\D y)^2})
$$
于是得到计算公式:
$$
u_{ij}^{k+1}=-u_{ij}^{k-1}+2u_{ij}^k+(\D t)^2\lambda(\frac{u_{i+1,j}^k-2u_{ij}^k+u_{i-1,j}^k}{\D x^2}+\frac{u_{i,j+1}^k-2u_{ij}^k+u_{i,j-1}^k}{\D y^2})
$$
初始条件的处理:${\partial u\over\partial t}|_{t=0}=0$,即
$$
u_{ij}^{-1}=u_{ij}^0
$$
边界条件的处理:$u|_{\text{boundary}}=0$,即
$$
u_{-1,j}^k =u_{N_x,j}^k= u_{i,-1}^k=u_{i,N_y}^k=0
$$
编写程序[T3.py](T3.py),获得解的动画[T3.mp4](./T3.mp4).
拖动视频进度条,寻找到振动从波峰到第一次回到水平位置的时间约为$t_0=0.2175$(以视频内显示的时间为准),代入计算得到:$\sqrt{5}\pi t_0=1.527\approx \frac \pi 2\approx 1.57$,计算值和理论值符合得很好
#### (3)与(4)
$$
\D t\le \frac{1}{\sqrt \lambda }\left(\frac{1}{\D x^2}+\frac{1}{\D y^2}\right)^{-1/2}
$$
取定$\D x=\D y=40$,$\lambda = 1$,得到稳定条件$\D t\le 0.0177$.上一问的模拟中$N_t=400$,时间0-3s,$\D t = \f 3 {400}=7.5\times 10^{-3}$,满足条件.

分别取$\D t=0.01,\D t=0.015,\D t=0.03$,模拟结果:分别如[T3_Dt=0.0100.mp4](./T3_Dt%3D0.0100.mp4),[T3_Dt=0.0150.mp4](./T3_Dt%3D0.0150.mp4),[T3_Dt=0.0300.mp4](./T3_Dt%3D0.0300.mp4)所示.其中当$0.03s$时,模拟过程中出现大量错误,可见数值稳定性条件的重要性.
```
Traceback (most recent call last):
  File "D:\anaconda\lib\site-packages\matplotlib\cbook\__init__.py", line 270, in process
    func(*args, **kwargs)
  File "D:\anaconda\lib\site-packages\matplotlib\backend_bases.py", line 3063, in mouse_move      
    s = self._mouse_event_to_message(event)      
  File "D:\anaconda\lib\site-packages\matplotlib\backend_bases.py", line 3055, in _mouse_event_to_message
    data_str = a.format_cursor_data(data).rstrip()
  File "D:\anaconda\lib\site-packages\matplotlib\image.py", line 1005, in format_cursor_data      
    self.colorbar.formatter.format_data_short(data)).strip()
  File "D:\anaconda\lib\site-packages\matplotlib\ticker.py", line 751, in format_data_short       
    - math.floor(math.log10(delta)))
ValueError: math domain error
```
上面是$\D x=\D y$的情况.对于$\D x\ne \D y$,验证如下:取$N_x=30,N_y=80$控制时间轴采样点数不变,总时长分别设为$3s,4s,5s$(其中5s下是不稳定的),得到的动画如视频所示.
[3s](T3_Dt=0.0075,Dx=0.033,Dy=0.013,stable.mp4)
[4s](T3_Dt=0.0100,Dx=0.033,Dy=0.013,stable.mp4)
[5s](T3_Dt=0.0125,Dx=0.033,Dy=0.013,unstable.mp4)

## 4.随机数的卡方检验

容易知道,落在各个区间的随机数频数值满足相同的概率分布.
进行多次如下实验:将$[0,1]$等分成$k$个区间,取$N$个在此区间均匀分布的随机数,随机数落在各个区间的数量分别为$n_{k'},k'=1,\cdots,k$,计算卡方值$\chi^2 \equiv \sum_{k'=1}^{k} \frac{(n_k'-\overline{n_k'})^2}{\overline{n_k'}}$,则该值服从参数为$k-1$的卡方分布:
$$
\chi^2\sim\chi^2(k-1)
$$
概率密度函数为:
$$
P_{v}(x)=\frac{1}{2^{v / 2} \Gamma(v / 2)}  t^{v / 2-1} e^{-t / 2} d t,\qquad v=k-1
$$
作程序[T4.py](T4.py),得到频率分布直方图和理论函数曲线如图所示:
\img["T4_1.png"]
可以看到,频率分布和理论曲线符合得很好.
#### (2)
随机选择一次实验的结果:$\chi^2=8.2$
以置信度$\a=0.05$检验:查表得到$\chi^2(9)|_{\a=0.05}=16.92\gt \chi^2$,可假定随机性正确
以置信度$\a=0.01$检验:查表得到$\chi^2(9)|_{\a=0.01}=21.67\gt \chi^2$,可假定随机性正确