# %% [markdown]
# 所求方程：
# $$
# \ddot{x} + \omega_0^2 x + Ax^2 +B x^3 + ... = 0
# $$
# 假设：
# $$
# x(t)=\sum_{n=-\infty}^{\infty} f_n e^{i n \omega t}\\
# \omega = \omega_0 + \omega_1 +\omega_2 +\cdots \\
# f_n=f_n^{(0)}+f_n^{(1)}+f_n^{(2)}+\cdots
# $$
# 用阶数命名AB等参数：
# $$
# \omega_0^2\to A_0\\
# A\to A_1 \\
# B\to A_2 \\
# C\to A_3
# $$

# %% [markdown]
# $$
# \omega_0^2 x <<Ax^2 <<Bx^3
# $$

# %%

import numpy,scipy,sympy
# from IPython.display import display, Markdown,Latex,display_latex
# from utils import *
import sympy
import itertools,functools
from sympy.abc import a,b,x,p,n

MAXLEN=10
iteration=8
def getidx(i):
    return i+MAXLEN

#---
import sys
sys.stdout=open("output.txt",'w')
def display(x):
    print(x)
def Latex(x):
    return x
def Markdown(x):
    return x
#---

#-------------符号定义----------------------
omegas=sympy.symbols(','.join(["omega_{}".format(i) for i in range(iteration)]),positive=True)
sym_f =[[sympy.symbols("f_{{{}}}^{{({})}}".format(i,j)) for i in range(-MAXLEN,MAXLEN+1)] for j in range(iteration)] #所有符号
f_sym=[sympy.symbols("f^{{({})}}".format(i)) for i in range(iteration)]
params=[omegas[0]**2]+[sympy.symbols("A_{{{}}}".format(i)) for i in range(1,iteration)] #Ax^2,Bx^3,...
res_f=[[0 for i in range(-MAXLEN,MAXLEN+1)] for j in range(iteration)] #存储结果
res_f[0][getidx(1)] = a       #初值
res_f[0][getidx(-1)] = a.conjugate()
res_conv_f =[ [[0 for i in range(-MAXLEN,MAXLEN+1)] for j in range(iteration)] for k in range(iteration)]
res_conv_f[0][0]=res_f[0][:]
# res_conv_f[i][k][getidx(n)] #i阶项,k次卷积，第n个值
#p 为阶数，计算时为了求出统一阶数，方便合并同阶项，用p为宗量构造多项式，
# 即 omega=omega_0 +omega_1 + ... ---->    omega_0 p^0 + omega_1 p^1 + ...
f_sym_tmp = sympy.Poly(reversed(f_sym),p)  #多项式，构造用
omega_tmp=sympy.Poly(reversed(omegas),p) #omega对应多项式

# 构建非线性项：
nonlin=0
for i in range(iteration):
    nonlin+=p**i * f_sym_tmp**(i+1) *params[i]
coeffs_tmp = nonlin-n**2 * f_sym_tmp*omega_tmp**2  #e^{in\omega t}的系数的表达式，对应n的每一项都应该为零
coeffs=[]
for i in range(iteration):
    coeffs.append(coeffs_tmp.coeff_monomial(p**i)) #coeffs[i]: p**i为 order i 的系数表达式，最高取到多少阶是准确的？iteration order！
# coeffs[3]

#-------------end-------------------------




#---------------------函数定义-----------------------------

def dconv(a,b,ret=None): #discrete convolution
    if ret==[]: 
        append,ref=True,True
    elif ret==None: #return the coeffsult. Do not save in ret
        append,ref=False,False
    else:
        append,ref=False,True
    tmp=[]
    for u in range(-MAXLEN,MAXLEN+1):
        s=0
        for n in range(max(u-MAXLEN,-MAXLEN),min(u+MAXLEN,MAXLEN)+1):
            s+= a[getidx(n)]*b[getidx(u-n)] 
        if ref:
            if append:
                ret.append(s)
            else:
                ret[getidx(u)]=s
        else:
            tmp.append(s)
    return tmp
def dconv_n(a,b,n): #计算离散卷积的第n项
    s=0
    assert n>=-MAXLEN and n<= MAXLEN
    for i in range(max(n-MAXLEN,-MAXLEN),min(n+MAXLEN,MAXLEN)+1):
        s+= a[getidx(i)]*b[getidx(n-i)]
    return s
def kPower_dconv_n(a,k,n):
    assert n>= -MAXLEN and n<= MAXLEN 
    cp=a[:]
    for i in range(k-1):
        cp=dconv(cp,a)
    return cp[n]
def gen_dconv(a):
    cp=a[:]
    while True:
        yield cp
        cp=dconv(cp,a)

def show(f,upperidx=None): #显示输出
    if not  upperidx: #不输出上标
        for i in range(MAXLEN+1):
            if f[getidx(i)]==0:
                continue
            if i==0:
                ltx=sympy.latex(f[getidx(0)])
                display(Latex(f"$f_0={ltx}$"))
            else:
                ltx=sympy.latex(f[getidx(i)])
                ltx2=sympy.latex(f[getidx(-i)])
                display(Latex(f"$f_{{{i}}}={ltx},\qquad f_{{{-i}}}={ltx2}$"))
    else:
        display(Markdown(f"# 第{upperidx}阶："))
        for i in range(MAXLEN+1):
            if f[getidx(i)]==0:
                continue
            if i==0:
                ltx=sympy.latex(f[getidx(0)])
                display(Latex(f"$f_0^{{{upperidx}}}={ltx}$"))
                # print(f"$f_0^{{{upperidx}}}={ltx}$")
            else:
                ltx=sympy.latex(f[getidx(i)])
                ltx2=sympy.latex(f[getidx(-i)])
                display(Latex(f"$f_{{{i}}}^{{{upperidx}}}={ltx},\qquad f_{{{-i}}}^{{{upperidx}}}={ltx2}$"))
                # print(f"$f_{{{i}}}^{{{upperidx}}}={ltx},\qquad f_{{{-i}}}^{{{upperidx}}}={ltx2}$")
def evaluate(expr):
    def addlist(*args):
        ret=[]
        
        for i in range(-MAXLEN,MAXLEN+1):
            s=0
            for ele in args:
                s+= ele[getidx(i)]
            ret.append(s)
            # ret.append(sum([ele[getidx(i)] for ele in args]))
        return ret 
    def eval(expr):
        if expr.func == sympy.Add: #应该为表达式树的最顶层，所有列表对应项加起来
            args=[]
            for i in range(len(expr.args)):
                args.append(eval(expr.args[i]))
            return addlist(*args)
        elif expr.func == sympy.Mul: #考虑(f^(k1))^a1 * (f^(k2))^a2 *... 项
            length= len(expr.args)
            is_f = [False for _ in range(length)]
            pre_calc=[None for _ in range(length)]
            for idx,arg in enumerate(expr.args):
                if arg.func == sympy.Symbol and arg in f_sym: #一次项
                    order=f_sym.index(arg)
                    is_f[idx]=True
                    pre_calc[idx]=res_conv_f[order][0][:]
                elif arg.func == sympy.Pow and arg.args[0] in f_sym: #高阶卷积项
                    order=f_sym.index(arg.args[0])
                    power=arg.args[1]
                    is_f[idx]=True
                    pre_calc[idx]=res_conv_f[order][power-1][:]
            to_do_conv_set=[]
            counter_set=[]
            for i in range(length):
                if is_f[i] :
                    to_do_conv_set.append(pre_calc[i])
                else:
                    counter_set.append(expr.args[i])
            ret = functools.reduce(lambda x,y:dconv(x,y),to_do_conv_set) #列表
            const=functools.reduce(lambda x,y:x*y,counter_set,1)
            return [const.subs(n,i) * ret[getidx(i)] for i in range(-MAXLEN,MAXLEN+1)]
        else:
            raise NotImplementedError("Not implemented.")
    return eval(expr)
def evaluate_omega(expr,symbol):
    k=expr.coeff(symbol).subs(n,1)
    b=(expr-expr.coeff(symbol)*symbol).subs(n,1)

    
    k_=evaluate(k)[getidx(1)] #(-n^2+1) \omega_0^2
    b_=evaluate(b)[getidx(1)]

    display(Markdown(f"$k\omega+b=0$，其中："))
    display(Latex(f"$$k={sympy.latex(k)}\\\\b={sympy.latex(b)}$$"))
    display(Latex(f"after evaluation:$$k={sympy.latex(k_)}\\\\b={sympy.latex(b_)}$$"))
    return sympy.simplify(-b_/k_)
def evaluate_f(expr,symbol):

    b=expr-expr.coeff(symbol)*symbol
    k_ = [(-i**2+1)*params[0] for i in range(-MAXLEN,MAXLEN+1)] #系数是固定的，每一阶都相同
    b_:list=evaluate(b)
    fn=[]
    for i in range(-MAXLEN,MAXLEN+1):
        idx=getidx(i)
        if abs(i)==1:
            fn.append(0)
        else:
            fn.append(sympy.expand(-b_[idx]/k_[idx]))
    return fn 
def evaluate_order(order):
    #计算前一项的所有卷积：
    display(Markdown(f"# 第{order}阶计算："))
    pre=order-1
    generator = gen_dconv(res_f[pre][:])
    for i in range(iteration):
        res_conv_f[pre][i] =next(generator)
    #计算omegas[order]的值 n=1
    display(Markdown("### $\omega$的计算："))
    omega=evaluate_omega(coeffs[order],omegas[order])
    display(Markdown(f"### $\omega_{order}={sympy.latex(omega)}$"))
    for i in range(iteration): #对于全部系数：
        coeffs[i] =coeffs[i].subs(omegas[order],omega) #代入omega值
    display(Markdown(f"### $f^{{({order})}}$的计算："))
    res_f[order]=evaluate_f(coeffs[order],f_sym[order])
    show(res_f[order])
    
    
#f0 的卷积除以f0不是降幂次！所以不能用solve直接除

#----------------end----------------------


# %%

# main:
for i in range(1,iteration):
    evaluate_order(i)

# %% [markdown]
# ### 第二阶计算(以下为手动计算，有错误,请勿参考）：
# $$
# \omega_2 = \frac{2A_1(f^{(1)}*f^{(0)})_{\pm 1}+A_2(f^{(0)})^3_{\pm 1}}{2\omega_0 f^{(0)}_{\pm 1}}=\frac{5A_1^2}{3\omega_0^3}(a^2+b^2)\\
# $$
# $$
# f^0_{\pm 1}=a\pm i b\\
# (f^0)^3_{\pm 1}=0\\
# (f^{1}*f^{0})_{\pm 1}=-\frac{5A_1(a^2+b^2)(a\pm i b)}{3\omega_0^2}
# $$
# 代入$\omega_2$的值：

# %%
# show(dconv(res_conv_f[1][0],res_conv_f[0][0]))
ret=dconv(dconv(res_f[0],res_f[0]),res_f[0])
show(ret)
#计算卷积的结果（a,b）：
conv=[[0 for _ in range(-MAXLEN,MAXLEN+1)] for i in range(iteration)]
dconv(res_f[0],res_f[0],conv[0]) #conv[i]为i+1次卷积得到的结果
for i in range(iteration-1):
    dconv(res_f[0],conv[i],conv[i+1])
#把最高阶的项表示出来，相当于求ax+b=0,这样当我们有了低阶项之后 可以一项一项地向上计算

fn_res=[]
# display(Markdown("# 求出各个f^i的表达式："))
for i in range(iteration):
    result = sympy.solve(coeffs[i],f_sym[i])[0]
    fn_res.append(result)
    # display(Latex(f"${sympy.latex(f_sym[i])}={sympy.latex(fn_res[i])}$"))
omega_res=[None]
#
# for i in range(1,iteration):# 第零阶不可计算
#     coeffs_1=coeffs[i].subs(n,1)# 代入n=1，计算omega的值
#     omega_res.append(sympy.solve(coeffs_1,omegas[i])[0])
    # display(Latex(f"${sympy.latex(omegas[i])}={sympy.latex(omega_res[i])}$")) 
#我们通过第一阶的计算，求出了omega_1=0，和$f_n^(1)$的值(n=+-1 undefined?)。用expr.subs来化简式子：

for i in range(1,iteration):
    coeffs_1 = coeffs[i].subs(n,1)
    omega_res.append(sympy.expand(sympy.solve(coeffs_1,omegas[i])[0]))
    display(Latex(f"${sympy.latex(omegas[i])}={sympy.latex(omega_res[i])}$")) 

# display(Markdown("# 代入$\omega_1=0$之后："))
res_f[1] = [params[1]/params[0]/(n**2-1) *dconv_n(res_f[0],res_f[0],n) if abs(n)!=1 else 0  for n in range(-MAXLEN,MAXLEN+1)]

for i in range(iteration):
    fn_res[i]=fn_res[i].subs(omegas[1],0) #化简
    # display(Latex(f"${sympy.latex( f_sym[i])}={sympy.latex(fn_res[i])}$"))
show(res_f[1])
for i in range(iteration):
    fn_res[i]=fn_res[i].subs(omegas[2],5*params[1]**2/3/omegas[0]**3*(a**2+b**2))
    # display(Latex(f"${sympy.latex(f_sym[i])}={sympy.latex(fn_res[i])}$"))


# %%
res_f[2]=[(2*params[1] *dconv_n(res_f[0],res_f[1],n)+ params[2] *conv[1][getidx(n)] -10*params[1]**2/(3*omegas[0]**2) *n**2*res_f[0][getidx(n)]*(a**2+b**2) )/((n**2-1)*params[0]) if abs(n)!=1 else 0 for n in range(-MAXLEN,MAXLEN+1)]
# show(coef[0],0)
# show(res_f[1],1)
show(res_f[2],2)

# %% [markdown]
# # 第三阶计算：
# 
# $$
# \omega_3=0
# $$

# %%
# show(conv[2]) +-1项为0
# show(dconv(res_f[1],res_f[1])) #+-1项为0
# show(dconv(res_f[0],res_f[2])) #0
# show(dconv(conv[0],res_f[1])) #0
# show(res_f[1]) #0
# show(res_f[0]) #not 0
for i in range(iteration):
    fn_res[i]=fn_res[i].subs(omegas[3],0)
    # display(Latex(f"${sympy.latex(f_sym[i])}={sympy.latex(fn_res[i])}$"))

# %%
res_f[3]=[(-(10*params[1]**2*res_f[1][getidx(n)]*n**2*(a**2+b**2))/(3*params[0])+2*params[1]*dconv_n(res_f[0],res_f[2],n)+params[1]*dconv_n(res_f[1],res_f[1],n)+3*params[2]*dconv_n(conv[0],res_f[1],n)+params[3]*conv[2][getidx(n)] )/(params[0]*(n**2-1)) if abs(n)!=1 else 0 for n in range(-MAXLEN,MAXLEN+1)]
show(res_f[3])

# %%
u,v=sympy.symbols("u v ")
w=sympy.expand(u*(v+1))
w=(u*(v+1))
sympy.srepr(w)
sympy.Pow(w,2)
expr = sympy.Add(u,u)
expr.func
sympy.core.numbers.Zero()
expr = 3*u **2 * v
expr
expr.func,expr.args
w=expr.args[2]
w.args[0].args


