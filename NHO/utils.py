import sympy
from sympy.abc import a,b,x,p,n
MAXLEN=10
iteration=7
def getidx(i):
    return i+MAXLEN
#-------------符号定义----------------------

omegas=sympy.symbols(','.join(["omega_{}".format(i) for i in range(iteration)]),positive=True)
sym_f =[[sympy.symbols("f_{{{}}}^{{({})}}".format(i,j)) for i in range(-MAXLEN,MAXLEN+1)] for j in range(iteration)] #所有符号
f_sym=[sympy.symbols("f^{{({})}}".format(i)) for i in range(iteration)]
params=[omegas[0]**2]+[sympy.symbols("A_{{{}}}".format(i)) for i in range(1,iteration)] #Ax^2,Bx^3,...
res_f=[[0 for i in range(-MAXLEN,MAXLEN+1)] for j in range(iteration)] #存储结果
res_f[0][getidx(1)] = a + 1j * b       #初值
res_f[0][getidx(-1)] = a - 1j * b
res_conv_f =[ [[0 for i in range(-MAXLEN,MAXLEN+1)] for j in range(iteration)] for k in range(iteration)]
res_conv_f[0][0]=res_f[0][:]
# res_conv_f[i][k][getidx(n)] #i阶项,k次卷积，第n个值
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
        print(cp)
        yield cp
        cp=dconv(cp,a)

def show(f,upperidx=None): #显示输出
    # display(Markdown())
    if not  upperidx: #不输出上标
        for i in range(MAXLEN+1):
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
            if i==0:
                ltx=sympy.latex(f[getidx(0)])
                display(Latex(f"$f_0^{{{upperidx}}}={ltx}$"))
            else:
                ltx=sympy.latex(f[getidx(i)])
                ltx2=sympy.latex(f[getidx(-i)])
                display(Latex(f"$f_{{{i}}}^{{{upperidx}}}={ltx},\qquad f_{{{-i}}}^{{{upperidx}}}={ltx2}$"))
def evaluate(expr,n):
    def eval(expr):
        func = expr.func 
        
        # if expr.func==sympy.Pow and expr.args[0] in f_sym:
        #     idx=f_sym.index(expr.args[0]) #第几阶项
        #     power = expr.args[1]
        #     return res_conv_f[idx][power-1][n]

        if expr.func == sympy.Mul: #考虑(f^(k1))^a1 * (f^(k2))^a2 *... 项
            is_f = [False for _ in range(len(expr.args))]
            pre_calc=[None for _ in range(len(expr.args))]
            for idx,arg in enumerate(expr.args):
                if arg.func == sympy.Symbol and arg in f_sym: #一次项
                    is_f[idx]=True
                    pre_calc[idx]=res_conv_f[idx][0][:]
                elif arg.func == sympy.Pow and arg.args[0] in f_sym:
                    power=expr.args[1]
                    is_f[idx]=True
                    pre_calc[idx]=res_conv_f[idx][power-1][:]
        if func == sympy.Symbol or func == sympy.Integer :
            return expr
        # return None
        args=list(map(evaluate,expr.args))
        return func(*args)
def evaluate_omega(expr,symbol):
    k=expr.coeff(symbol)
    b=expr-expr.coeff(symbol)*symbol
    evaluate(expr=k)
def evaluate_order(order):
    #计算前一项的所有卷积：
    pre=order-1
    generator = gen_dconv(res_f[pre][:])
    for i in range(iteration):
        res_conv_f[pre][i] =next(generator)
    evaluate_omega(coeffs[order],omegas[order])
    # Walking the Tree:
    # def pre(expr):
    #     display(Latex(f"${sympy.latex(expr)}$"))
    #     print(len(expr.args),expr.func)
    #     for arg in expr.args:
    #         pre(arg)

#f0 的卷积除以f0不是降幂次！所以不能用solve直接除

#----------------end----------------------
