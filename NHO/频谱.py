# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
MAXLEN=6 #列表项从n=-10到n=+10
from cycler import cycler


# 全局修改rcParams["axes.prop_cycle"]，对所有子图生效
default_cycler = (cycler(color=['r', 'g', 'b', 'y']) )
plt.rc('axes', prop_cycle=default_cycler)

def getidx(i):
    return i+MAXLEN

def draw(serial,alpha=1):
    '''
    serial是离散的频率谱
    '''
    ax=plt.gca()
    x=[alpha*i for i in range(-MAXLEN,MAXLEN+1)]
    # print(x,serial)
    # plt.set_cmap('Blues')
    
    
    plt.stem(x,np.array(serial)+np.array([0.01*i for i in range(-MAXLEN,MAXLEN+1)]),linefmt='C1-',markerfmt='C1o')
    plt.stem(x,serial)
    x_major_locator=plt.MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    plt.show()
# 一个使用示例：
s=[0.0 for i in range(-MAXLEN,MAXLEN+1)] #初始化一个频谱

s[getidx(1)]=s[getidx(-1)]=1 # 将n=+-1的值设为1
draw(s,1.1) #化s的频谱图
