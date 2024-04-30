#@title mybets { form-width: "100px" }

#from pyomo.core.expr import current as EXPR
from __future__ import division
import pyomo.core.expr as EXPR
import pyomo.environ as pyo
from pyomo.environ import *

def mybets(m):

    c=[[],[]]
    Q=[]
    zval=-m.beta
    for i in range(m.s.value): # initialise internal stability functions
        Q.append(m.b[i]*zval)

    for J in range(1,m.s.value): # determine A^{J=1..s-1}=>beta_{J+1} and store in inew
    #for J in range(1,int(m.order)): # determine A^{J=1..s-1}=>beta_{J+1} and store in inew
    #for J in range(0,1): # determine J+1 and store in inew
        ibet=J+1
        iold=(J-1)%2
        inew=J%2
        #print(iold,inew)
        for i in range(J):  #dead rows
            c[inew].append([None])

        #Populate c[inew]
        for i in range(J,m.s.value):  #live rows in A^J
            r1=[]
            for j in range(i-J+1): # live cols in A^J
                if J>1:
                    c1=0
                    for k in range(j+1,i-J+2): # non-zero vals: c[iold] in row i has non zero up to c[iold][i][i-J+1], a in col j has non zero from a[j+1,j]
                        #print(i,j,k)
                        c1+=c[iold][i][k]*m.a[k,j]
                else:
                    c1=m.a[i,j]   # initialise to A
                r1.append(c1)
            c[inew].append(r1)

        #Reset c[iold]
        c[iold]=[]

        #Use c[inew] to determine beta expression: first sum over row and then apply beta
        #csos23  ARE THESE MAD??
        if ibet>m.order:
            y1=0 # bet constraint expression
            for i in range(J,m.s.value): # live rows in A^J:   J to s-1 live rows (s-J live)
                rowsum=0
                for j in range(i-J+1):  # live cols in row i  of A^J (s-J-j live entries)
                    rowsum+=c[inew][i][j]
                y1+=m.b[i]*rowsum
            y1-=m.bet[ibet]

            m.cons.add(y1 == 0 )

            #print('**beta**',ibet,m.bet[ibet],EXPR.expression_to_string(y1),'==0\n')
            print('**beta**',ibet,m.bet[ibet])

        #Use c[inew] to determine internal stability functions evaluated at zval
        zpow=zval**(J+1)
        for j in range(m.s.value-J): # live cols in A^J:   0 to s-J-1 live cols (s-J live)
            colsum=0
            for i in range(J+j,m.s.value):  # live rows in col j  of A^J (s-J-j live entries)
                colsum+=m.b[i]*c[inew][i][j]
            Q[j]+=colsum*zpow #accumulate stability polynomials


    #csos23
    norm=0.0
    L=m.s.value
    #L=2
    for i in range(m.s.value): # initialise internal stability objective function
    #g=6
    #for i in range(g,g+1): # initialise internal stability objective function
        q1=abs(Q[i])
        #q1=Q[i]
        norm+=q1**L
    ##csos23
    ob=norm**(1/L)
    #ob=norm
    #print('Q',i,EXPR.expression_to_string(Q[i]),'\n')


    ##ALTERNATE OBJECTIVE FUNCTION
    norm=0.0
    L=m.s.value
    #L=2
    for i in range(1,m.s.value):
        for j in range(i):
            a1=abs(m.a[i,j])
            norm+=a1**L
    for j in m.cols:
        b1=abs(m.b[j])
        norm+=b1**L
    #csos
    ob=norm**(1/L)


    #ASSIGNING ob to objective function
    m.obj = pyo.Objective(expr=ob)

#@title mycons { form-width: "100px" }

# order conditions

import pyomo.environ as pyo
def mycons(m):

    c=[0]
    c.append(sum(m.a[1,i] for i in range(1)))
    c.append(sum(m.a[2,i] for i in range(2)))
    c.append(sum(m.a[3,i] for i in range(3)))
    c.append(sum(m.a[4,i] for i in range(4)))
    c.append(sum(m.a[5,i] for i in range(5)))
    c.append(sum(m.a[6,i] for i in range(6)))
    c.append(sum(m.a[7,i] for i in range(7)))
    c.append(sum(m.a[8,i] for i in range(8)))
    c.append(sum(m.a[9,i] for i in range(9)))
    c.append(sum(m.a[10,i] for i in range(10)))
    c.append(sum(m.a[11,i] for i in range(11)))
    c.append(sum(m.a[12,i] for i in range(12)))
    c.append(sum(m.a[13,i] for i in range(13)))
    c.append(sum(m.a[14,i] for i in range(14)))
    a=[0.03493615984916687, 0.03232833370566368, 1.0940261141456631e-08, 0.04490070790052414, 0.03454554080963135, 0.04777379333972931, 0.11284931749105453, 0.051660437136888504, 0.04164694622159004, 0.049222785979509354, 0.05507682263851166, 0.05214410275220871, 0.0430271215736866, 1.118902659413834e-08, 4.978309764425148e-09, 0.05601559206843376, 0.06165182217955589, 0.06048295646905899, 0.05834311619400978, 0.04131209850311279, 0.04766211658716202, 0.05363506078720093, 0.05120457708835602, 0.05579254776239395, 0.05486460030078888, 0.05182594060897827, 0.047686364501714706, 0.041190002113580704, 0.04811837524175644, 0.05222948640584946, 0.056764136999845505, 0.05625079572200775, 0.05472683906555176, 0.05594264343380928, 0.04189429432153702, 0.04112095758318901, 0.05002988874912262, 0.055259622633457184, 0.05975879356265068, 0.05313216149806976, 0.06190238147974014, 0.06198401376605034, 0.06360447406768799, 0.0465279147028923, 0.04958430677652359, 0.04079129546880722, 0.04715939238667488, 0.04877427965402603, 0.04635829105973244, 0.053854912519454956, 0.06270197778940201, 0.06563377380371094, 0.05587083473801613, 0.05359280854463577, 0.05121840909123421, 0.03374172002077103, 0.037938766181468964, 0.04218684136867523, 0.04757606238126755, 0.05292386934161186, 0.061214957386255264, 0.06253129988908768, 0.05774383246898651, 0.05909852311015129, 0.05504361912608147, 0.048739220947027206, 0.029355555772781372, 0.03562483564019203, 0.04547969624400139, 0.04715670645236969, 0.048295002430677414, 0.05496576800942421, 0.061271604150533676, 0.06462185829877853, 0.06177794933319092, 0.07091038674116135, 0.050187304615974426, 0.0490896999835968, 0.03606441617012024, 0.03483543545007706, 0.03423730656504631, 0.03898352012038231, 0.04096708446741104, 0.04637855291366577, 0.05352979898452759, 0.0555088147521019, 0.06145523488521576, 0.06719944626092911, 0.06396143138408661, 0.05609923228621483, 0.05008354038000107, 0.03547520563006401, 0.031185569241642952, 0.030588814988732338, 0.02862359955906868, 0.04091302305459976, 0.039865173399448395, 0.05024731159210205, 0.05380551517009735, 0.05950777977705002, 0.07241976261138916, 0.07978709042072296, 0.07425735145807266, 0.07512682676315308, 0.06416870653629303]
    b=[2.0855532856245418e-09, 3.69855657211815e-09, 0.16303128004074097, 1.0262908745062305e-08, 1.0272779071840432e-08, 0.15383031964302063, 0.14535488188266754, 1.3887297978243396e-09, 7.659721745767456e-09, 2.7695363780111393e-09, 2.9388642630578943e-09, 5.337066788513312e-09, 1.0660181182231554e-08, 0.23244597017765045, 0.27846816182136536]
    

    m.a[1,0] = a[0]
    m.a[2,0] = a[1]
    m.a[2,1] = a[2]
    m.a[3,0] = a[3]
    m.a[3,1] = a[4]
    m.a[3,2] = a[5]
    m.a[4,0] = a[6]
    m.a[4,1] = a[7]
    m.a[4,2] = a[8]
    m.a[4,3] = a[9]
    m.a[5,0] = a[10]
    m.a[5,1] = a[11]
    m.a[5,2] = a[12]
    m.a[5,3] = a[13]
    m.a[5,4] = a[14]
    m.a[6,0] = a[15]
    m.a[6,1] = a[16]
    m.a[6,2] = a[17]
    m.a[6,3] = a[18]
    m.a[6,4] = a[19]
    m.a[6,5] = a[20]
    m.a[7,0] = a[21]
    m.a[7,1] = a[22]
    m.a[7,2] = a[23]
    m.a[7,3] = a[24]
    m.a[7,4] = a[25]
    m.a[7,5] = a[26]
    m.a[7,6] = a[27]
    m.a[8,0] = a[28]
    m.a[8,1] =  a[29]
    m.a[8,2]  = a[30]
    m.a[8,3] = a[31]
    m.a[8,4] = a[32]
    m.a[8,5] = a[33]
    m.a[8,6] = a[34]
    m.a[8,7] = a[35]
    m.a[9,0] =  a[36]
    m.a[9,1] = a[37]
    m.a[9,2] = a[38]
    m.a[9,3] = a[39]
    m.a[9,4] = a[40]
    m.a[9,5] =  a[41]
    m.a[9,6] = a[42]
    m.a[9,7] = a[43]
    m.a[9,8] = a[44]
    m.a[10,0] = a[45]
    m.a[10,1] = a[46]
    m.a[10,2] = a[47]
    m.a[10,3] = a[48]
    m.a[10,4] = a[49]
    m.a[10,5] = a[50]
    m.a[10,6] = a[51]
    m.a[10,7] = a[52]
    m.a[10,8] = a[53]
    m.a[10,9] = a[54]
    m.a[11,0] = a[55]
    m.a[11,1] = a[56]
    m.a[11,2] = a[57]
    m.a[11,3] = a[58]
    m.a[11,4] = a[59]
    m.a[11,5] = a[60]
    m.a[11,6] = a[61]
    m.a[11,7] = a[62]
    m.a[11,8] = a[63]
    m.a[11,9] = a[64]
    m.a[11,10]= a[65]
    m.a[12,0] = a[66]
    m.a[12,1] = a[67]
    m.a[12,2] = a[68]
    m.a[12,3] = a[69]
    m.a[12,4] = a[70]
    m.a[12,5] = a[71]
    m.a[12,6] = a[72]
    m.a[12,7] = a[73]
    m.a[12,8] = a[74]
    m.a[12,9] = a[75]
    m.a[12,10]= a[76]
    m.a[12,11]= a[77]
    m.a[13,0] = a[78]
    m.a[13,1] = a[79]
    m.a[13,2] = a[80]
    m.a[13,3] = a[81]
    m.a[13,4] = a[82]
    m.a[13,5] = a[83]
    m.a[13,6] = a[84]
    m.a[13,7] = a[85]
    m.a[13,8] = a[86]
    m.a[13,9] = a[87]
    m.a[13,10] =a[88]
    m.a[13,11] =a[89]
    m.a[13,12] =a[90]
    m.a[14,0] = a[91]
    m.a[14,1] = a[92]
    m.a[14,2] = a[93]
    m.a[14,3] = a[94]
    m.a[14,4] = a[95]
    m.a[14,5] = a[96]
    m.a[14,6] = a[97]
    m.a[14,7] = a[98]
    m.a[14,8] = a[99]
    m.a[14,9] = a[100]
    m.a[14,10] = a[101]
    m.a[14,11] = a[102]
    m.a[14,12] = a[103]
    m.a[14,13] = a[104]


 


    

     

    m.b[0] = b[0]
    m.b[1] = b[1]
    m.b[2] = b[2]
    m.b[3] = b[3]
    m.b[4] = b[4]
    m.b[5] = b[5]
    m.b[6] = b[6]
    m.b[7] = b[7]
    m.b[8] = b[8]
    m.b[9] = b[9]
    m.b[10] = b[10]
    m.b[11] = b[11]
    m.b[12] = b[12]
    m.b[13] = b[13]
    m.b[14] = b[14]

    y1 = +m.b[0]+m.b[1]+m.b[2]+m.b[3]+m.b[4]+m.b[5]+m.b[6]+m.b[7]+m.b[8]+m.b[9]+m.b[10]+m.b[11]+m.b[12]+m.b[13]+m.b[14] - (1/1)
    m.cons.add(y1 == 0 )

    y1 = +m.b[1]*c[1]+m.b[2]*c[2]+m.b[3]*c[3]+m.b[4]*c[4]+m.b[5]*c[5]+m.b[6]*c[6]+m.b[7]*c[7]+m.b[8]*c[8]+m.b[9]*c[9]+m.b[10]*c[10]+m.b[11]*c[11]+m.b[12]*c[12]+m.b[13]*c[13]+m.b[14]*c[14] - (1/2)
    m.cons.add(y1 == 0 )

    y1 = +m.b[1]*c[1]**2+m.b[2]*c[2]**2+m.b[3]*c[3]**2+m.b[4]*c[4]**2+m.b[5]*c[5]**2+m.b[6]*c[6]**2+m.b[7]*c[7]**2+m.b[8]*c[8]**2+m.b[9]*c[9]**2+m.b[10]*c[10]**2+m.b[11]*c[11]**2+m.b[12]*c[12]**2+m.b[13]*c[13]**2+m.b[14]*c[14]**2 - (1/3)
    m.cons.add(y1 == 0 )

    y1 = +m.b[2]*m.a[2,1]*c[1]+m.b[3]*m.a[3,1]*c[1]+m.b[3]*m.a[3,2]*c[2]+m.b[4]*m.a[4,1]*c[1]+m.b[4]*m.a[4,2]*c[2]+m.b[4]*m.a[4,3]*c[3]+m.b[5]*m.a[5,1]*c[1]+m.b[5]*m.a[5,2]*c[2]+m.b[5]*m.a[5,3]*c[3]+m.b[5]*m.a[5,4]*c[4]+m.b[6]*m.a[6,1]*c[1]+m.b[6]*m.a[6,2]*c[2]+m.b[6]*m.a[6,3]*c[3]+m.b[6]*m.a[6,4]*c[4]+m.b[6]*m.a[6,5]*c[5]+m.b[7]*m.a[7,1]*c[1]+m.b[7]*m.a[7,2]*c[2]+m.b[7]*m.a[7,3]*c[3]+m.b[7]*m.a[7,4]*c[4]+m.b[7]*m.a[7,5]*c[5]+m.b[7]*m.a[7,6]*c[6]+m.b[8]*m.a[8,1]*c[1]+m.b[8]*m.a[8,2]*c[2]+m.b[8]*m.a[8,3]*c[3]+m.b[8]*m.a[8,4]*c[4]+m.b[8]*m.a[8,5]*c[5]+m.b[8]*m.a[8,6]*c[6]+m.b[8]*m.a[8,7]*c[7]+m.b[9]*m.a[9,1]*c[1]+m.b[9]*m.a[9,2]*c[2]+m.b[9]*m.a[9,3]*c[3]+m.b[9]*m.a[9,4]*c[4]+m.b[9]*m.a[9,5]*c[5]+m.b[9]*m.a[9,6]*c[6]+m.b[9]*m.a[9,7]*c[7]+m.b[9]*m.a[9,8]*c[8]+m.b[10]*m.a[10,1]*c[1]+m.b[10]*m.a[10,2]*c[2]+m.b[10]*m.a[10,3]*c[3]+m.b[10]*m.a[10,4]*c[4]+m.b[10]*m.a[10,5]*c[5]+m.b[10]*m.a[10,6]*c[6]+m.b[10]*m.a[10,7]*c[7]+m.b[10]*m.a[10,8]*c[8]+m.b[10]*m.a[10,9]*c[9]+m.b[11]*m.a[11,1]*c[1]+m.b[11]*m.a[11,2]*c[2]+m.b[11]*m.a[11,3]*c[3]+m.b[11]*m.a[11,4]*c[4]+m.b[11]*m.a[11,5]*c[5]+m.b[11]*m.a[11,6]*c[6]+m.b[11]*m.a[11,7]*c[7]+m.b[11]*m.a[11,8]*c[8]+m.b[11]*m.a[11,9]*c[9]+m.b[11]*m.a[11,10]*c[10]+m.b[12]*m.a[12,1]*c[1]+m.b[12]*m.a[12,2]*c[2]+m.b[12]*m.a[12,3]*c[3]+m.b[12]*m.a[12,4]*c[4]+m.b[12]*m.a[12,5]*c[5]+m.b[12]*m.a[12,6]*c[6]+m.b[12]*m.a[12,7]*c[7]+m.b[12]*m.a[12,8]*c[8]+m.b[12]*m.a[12,9]*c[9]+m.b[12]*m.a[12,10]*c[10]+m.b[12]*m.a[12,11]*c[11]+m.b[13]*m.a[13,1]*c[1]+m.b[13]*m.a[13,2]*c[2]+m.b[13]*m.a[13,3]*c[3]+m.b[13]*m.a[13,4]*c[4]+m.b[13]*m.a[13,5]*c[5]+m.b[13]*m.a[13,6]*c[6]+m.b[13]*m.a[13,7]*c[7]+m.b[13]*m.a[13,8]*c[8]+m.b[13]*m.a[13,9]*c[9]+m.b[13]*m.a[13,10]*c[10]+m.b[13]*m.a[13,11]*c[11]+m.b[13]*m.a[13,12]*c[12]+m.b[14]*m.a[14,1]*c[1]+m.b[14]*m.a[14,2]*c[2]+m.b[14]*m.a[14,3]*c[3]+m.b[14]*m.a[14,4]*c[4]+m.b[14]*m.a[14,5]*c[5]+m.b[14]*m.a[14,6]*c[6]+m.b[14]*m.a[14,7]*c[7]+m.b[14]*m.a[14,8]*c[8]+m.b[14]*m.a[14,9]*c[9]+m.b[14]*m.a[14,10]*c[10]+m.b[14]*m.a[14,11]*c[11]+m.b[14]*m.a[14,12]*c[12]+m.b[14]*m.a[14,13]*c[13] - (1/6)
    m.cons.add(y1 == 0 )



#@title beta values { form-width: "100px" }

# used for linear order conditions plus linear stability polynomial matching conditions

import numpy as np

data=np.array([
    72.362092230726359621444062630453895436513878479204372239744118914339008361957274894021123680041692242,
1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
0.50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
0.166666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667,
0.0330278406028757527979332298466336725484247797530949148541780288765549097681738463834524909538164545,
0.00407886754615175430254075890828995614461800841500213845583064974804470591686586532729380258593128180,
0.000333527373172737242603705467098764792409315403431559178662909610705332831105918682590799425317677283,
0.00001884387445457719649926991124773774639255686232669238908076544455395077431328132833751461219355974457,
7.55565600673228187844472799772168548435456673497418681132813388498239377614064205083716425567987086e-7,
2.180674889659337806106474225333623253164152989661516940939984945939720123180192751320554488732457046e-8,
4.54329575207387209531051040170803945739047297526622495153899072480429169194765507413726877684905152e-10,
6.76964099198340395235741216662494203744813306696814739898497900746954960650250467206606180555417237e-12,
7.03301065055001143525874293209920770968902766712091555585888491039951221050037050929236960249104521e-14,
4.83760685817544104974918008895344303724730995701958486014397001982907347822621723661931383126126731e-16,
1.979536791110686622489750097380272593062776391994465218701741931323003217460574992708738821233104875e-18,
3.64746555383158031489459760645481664971731724579909139392166493612827917797499294769223945733650777e-21
              ])



#@title solve { form-width: "100px" }
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverFactory,SolverStatus,TerminationCondition
#from pyomo.core.expr import current as EXPR
import pyomo.core.expr as EXPR
from pyomo.core.expr import value
import time
from random import uniform,seed,randrange
import sys

reloadmodel=False

#data = np.loadtxt('beta.dat')
betaval = data[0]
betdic = dict(enumerate(data[1:]))

minbet=min(data)
#csos23 bettol=1e-9*minbet
#csos23
bettol=1e-6 ## Could push to 15?


m = pyo.ConcreteModel()


##NB EDIT THIS BY HAND
m.order = pyo.Param(initialize=3)
m.macc = pyo.Param(initialize=5)
print('Order:',m.order.value,'Acc:',m.macc.value)



#csos23 m.s = pyo.Param(initialize=len(betdic)-1)
m.s = pyo.Param(initialize=m.order.value*m.macc.value)


initval= 1/m.s.value
lneg=0 # Set to 0 to disable negative coefficients in bounding and initialisation: 0 appropriate for SSP
lim=1
skip=1

print('lim:',lim,'skip:',skip,'initval:',initval,'lneg:',lneg)

def setlsrk(m):
    for j in range(m.s.value-skip): #cols
        for i in range(j+skip+1,m.s.value):
            m.cons.add(m.a[i,j] == m.b[j] )


def azeros(m):
    for i in m.rows:
        for j in range(i,m.s.value): #cols
            m.cons.add(m.a[i,j] == 0 )


m.rows= range(m.s.value)
m.cols= range(m.s.value)

m.a= pyo.Var( m.rows, m.cols, bounds=(-lneg*lim,lim))
m.b= pyo.Var( m.rows, bounds=(-lneg*lim,lim))
m.beta=pyo.Param(initialize=betaval)  #beta_R
m.bet = pyo.Param(range(m.s.value+1), initialize=betdic) #beta_j

print('beta_R:',m.beta.value)

m.cons = pyo.ConstraintList()

mycons(m)
#csos23
setlsrk(m)
#csos23?
#azeros(m) #Nones faster
mybets(m)

#m.pprint()


seedval = randrange(sys.maxsize)
#seedval=2842383511757470787
seed(seedval)
print('seedval:',seedval)



### solver = SolverFactory('ipopt',executable='/content/Ipopt-3.12.13/build/bin/ipopt')
solver = SolverFactory('ipopt',executable='ipopt')
#solver.options['max_iter'] = 10000
solver.options['nlp_scaling_method'] = 'none'
solver.options['linear_solver'] = 'ma27'
solver.options['honor_original_bounds']='yes'

solver.options['bound_relax_factor']=bettol ## default 1e-8   https://kiwi.oden.utexas.edu/papers/Multi-output-multi-fidelity-MLBLUE-Croci-Wright-Willcox.pdf

solver.options['tol'] = bettol # default 1e-8
solver.options['dual_inf_tol']=1e8*bettol #default 1, Larger more accurate
solver.options['constr_viol_tol']=1e4*bettol #default 1e-4
solver.options['compl_inf_tol']=1e4*bettol #default 1e-4

solver.options['acceptable_tol'] = bettol # default 1e-6
solver.options['acceptable_dual_inf_tol']=1e8*bettol #default 1e10
solver.options['acceptable_constr_viol_tol']=1e4*bettol #default 1e-2
solver.options['acceptable_compl_inf_tol']=1e4*bettol #default 1e-2



niter=1

#solver = SolverFactory('multistart')
#niter=1
#nitermulti=10

#solver = SolverFactory('couenne')
#niter=1


if reloadmodel:
    with open('m.pkl', mode='rb') as file: #reload model
        m = cloudpickle.load(file)
        print('Reload')
else:
    objvalbest=1e10
    for iter in range(niter):

        start = time.time()

        for i in range(1,m.s.value):
            for j in range(i):
                #csos23
                m.a[i,j].value=uniform(-lneg*initval,initval)
                #m.a[i,j].value=initval
        for j in m.cols:
            #csos23
            m.b[i].value=uniform(-lneg*initval,initval)
            #m.b[i].value=initval   ## sum b=1 !



        try:

            results=solver.solve(m,tee=True)
            #results=solver.solve(m,tee=False)

            #results=solver.solve(m,solver='ipopt',strategy='rand_distributed',iterations=nitermulti)
            if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
                end = time.time()
                objval=m.obj()
                print(iter,'Time:',end - start,'Objective:',objval)
                if objval < objvalbest:
                    objvalbest=objval
                    mbest=m.clone()
                    resultsbest=results
            elif (results.solver.termination_condition == TerminationCondition.infeasible):
                print(iter,'Infeasible.')
            else:
                print(iter,'Failed. Solver Status:',  result.solver.status)
        except:
            print(iter,'Ipopt error, moving on.')

        #m.display()
        #results.write()


    m=mbest.clone()  # revert to best instance
    results=resultsbest


    import cloudpickle  # write model to file
    with open('m.pkl', mode='wb') as file:
        cloudpickle.dump(m, file)


objval=m.obj()
print('Final objective value:',objval)

#@title stability polynomials { form-width: "100px" }


a = []
for i in range(1,m.s.value):
    for j in range(i):
              
            c = ('%5.20f'%m.a[i,j].value)
            a.append(float(c))

#print("\n")
b = []
#print("B")
for j in m.cols:
    c = ('%5.20f'%m.b[j].value)
    b.append(float(c))

print(a)
print(b)

	
#@title mybets { form-width: "100px" }

#from pyomo.core.expr import current as EXPR
from __future__ import division
import pyomo.core.expr as EXPR
import pyomo.environ as pyo
from pyomo.environ import *

def mybets(m):

    c=[[],[]]
    Q=[]
    zval=-m.beta
    for i in range(m.s.value): # initialise internal stability functions
        Q.append(m.b[i]*zval)

    for J in range(1,m.s.value): # determine A^{J=1..s-1}=>beta_{J+1} and store in inew
    #for J in range(1,int(m.order)): # determine A^{J=1..s-1}=>beta_{J+1} and store in inew
    #for J in range(0,1): # determine J+1 and store in inew
        ibet=J+1
        iold=(J-1)%2
        inew=J%2
        #print(iold,inew)
        for i in range(J):  #dead rows
            c[inew].append([None])

        #Populate c[inew]
        for i in range(J,m.s.value):  #live rows in A^J
            r1=[]
            for j in range(i-J+1): # live cols in A^J
                if J>1:
                    c1=0
                    for k in range(j+1,i-J+2): # non-zero vals: c[iold] in row i has non zero up to c[iold][i][i-J+1], a in col j has non zero from a[j+1,j]
                        #print(i,j,k)
                        c1+=c[iold][i][k]*m.a[k,j]
                else:
                    c1=m.a[i,j]   # initialise to A
                r1.append(c1)
            c[inew].append(r1)

        #Reset c[iold]
        c[iold]=[]

        #Use c[inew] to determine beta expression: first sum over row and then apply beta
        #csos23  ARE THESE MAD??
        if ibet>m.order:
            y1=0 # bet constraint expression
            for i in range(J,m.s.value): # live rows in A^J:   J to s-1 live rows (s-J live)
                rowsum=0
                for j in range(i-J+1):  # live cols in row i  of A^J (s-J-j live entries)
                    rowsum+=c[inew][i][j]
                y1+=m.b[i]*rowsum
            y1-=m.bet[ibet]

            m.cons.add(y1 == 0 )

            #print('**beta**',ibet,m.bet[ibet],EXPR.expression_to_string(y1),'==0\n')
            print('**beta**',ibet,m.bet[ibet])

        #Use c[inew] to determine internal stability functions evaluated at zval
        zpow=zval**(J+1)
        for j in range(m.s.value-J): # live cols in A^J:   0 to s-J-1 live cols (s-J live)
            colsum=0
            for i in range(J+j,m.s.value):  # live rows in col j  of A^J (s-J-j live entries)
                colsum+=m.b[i]*c[inew][i][j]
            Q[j]+=colsum*zpow #accumulate stability polynomials


    #csos23
    norm=0.0
    L=m.s.value
    #L=2
    for i in range(m.s.value): # initialise internal stability objective function
    #g=6
    #for i in range(g,g+1): # initialise internal stability objective function
        q1=abs(Q[i])
        #q1=Q[i]
        norm+=q1**L
    ##csos23
    ob=norm**(1/L)
    #ob=norm
    #print('Q',i,EXPR.expression_to_string(Q[i]),'\n')


    ##ALTERNATE OBJECTIVE FUNCTION
    norm=0.0
    L=m.s.value
    #L=2
    for i in range(1,m.s.value):
        for j in range(i):
            a1=abs(m.a[i,j])
            norm+=a1**L
    for j in m.cols:
        b1=abs(m.b[j])
        norm+=b1**L
    #csos
    ob=norm**(1/L)


    #ASSIGNING ob to objective function
    m.obj = pyo.Objective(expr=ob)

#@title mycons { form-width: "100px" }

# order conditions

import pyomo.environ as pyo
def mycons(m):

    c=[0]
    c.append(sum(m.a[1,i] for i in range(1)))
    c.append(sum(m.a[2,i] for i in range(2)))
    c.append(sum(m.a[3,i] for i in range(3)))
    c.append(sum(m.a[4,i] for i in range(4)))
    c.append(sum(m.a[5,i] for i in range(5)))
    c.append(sum(m.a[6,i] for i in range(6)))
    c.append(sum(m.a[7,i] for i in range(7)))
    c.append(sum(m.a[8,i] for i in range(8)))
    c.append(sum(m.a[9,i] for i in range(9)))
    c.append(sum(m.a[10,i] for i in range(10)))
    c.append(sum(m.a[11,i] for i in range(11)))
    c.append(sum(m.a[12,i] for i in range(12)))
    c.append(sum(m.a[13,i] for i in range(13)))
    c.append(sum(m.a[14,i] for i in range(14)))
    a=[0.03493615984916687, 0.03232833370566368, 1.0940261141456631e-08, 0.04490070790052414, 0.03454554080963135, 0.04777379333972931, 0.11284931749105453, 0.051660437136888504, 0.04164694622159004, 0.049222785979509354, 0.05507682263851166, 0.05214410275220871, 0.0430271215736866, 1.118902659413834e-08, 4.978309764425148e-09, 0.05601559206843376, 0.06165182217955589, 0.06048295646905899, 0.05834311619400978, 0.04131209850311279, 0.04766211658716202, 0.05363506078720093, 0.05120457708835602, 0.05579254776239395, 0.05486460030078888, 0.05182594060897827, 0.047686364501714706, 0.041190002113580704, 0.04811837524175644, 0.05222948640584946, 0.056764136999845505, 0.05625079572200775, 0.05472683906555176, 0.05594264343380928, 0.04189429432153702, 0.04112095758318901, 0.05002988874912262, 0.055259622633457184, 0.05975879356265068, 0.05313216149806976, 0.06190238147974014, 0.06198401376605034, 0.06360447406768799, 0.0465279147028923, 0.04958430677652359, 0.04079129546880722, 0.04715939238667488, 0.04877427965402603, 0.04635829105973244, 0.053854912519454956, 0.06270197778940201, 0.06563377380371094, 0.05587083473801613, 0.05359280854463577, 0.05121840909123421, 0.03374172002077103, 0.037938766181468964, 0.04218684136867523, 0.04757606238126755, 0.05292386934161186, 0.061214957386255264, 0.06253129988908768, 0.05774383246898651, 0.05909852311015129, 0.05504361912608147, 0.048739220947027206, 0.029355555772781372, 0.03562483564019203, 0.04547969624400139, 0.04715670645236969, 0.048295002430677414, 0.05496576800942421, 0.061271604150533676, 0.06462185829877853, 0.06177794933319092, 0.07091038674116135, 0.050187304615974426, 0.0490896999835968, 0.03606441617012024, 0.03483543545007706, 0.03423730656504631, 0.03898352012038231, 0.04096708446741104, 0.04637855291366577, 0.05352979898452759, 0.0555088147521019, 0.06145523488521576, 0.06719944626092911, 0.06396143138408661, 0.05609923228621483, 0.05008354038000107, 0.03547520563006401, 0.031185569241642952, 0.030588814988732338, 0.02862359955906868, 0.04091302305459976, 0.039865173399448395, 0.05024731159210205, 0.05380551517009735, 0.05950777977705002, 0.07241976261138916, 0.07978709042072296, 0.07425735145807266, 0.07512682676315308, 0.06416870653629303]
    b=[2.0855532856245418e-09, 3.69855657211815e-09, 0.16303128004074097, 1.0262908745062305e-08, 1.0272779071840432e-08, 0.15383031964302063, 0.14535488188266754, 1.3887297978243396e-09, 7.659721745767456e-09, 2.7695363780111393e-09, 2.9388642630578943e-09, 5.337066788513312e-09, 1.0660181182231554e-08, 0.23244597017765045, 0.27846816182136536]
    

    m.a[1,0] = a[0]
    m.a[2,0] = a[1]
    m.a[2,1] = a[2]
    m.a[3,0] = a[3]
    m.a[3,1] = a[4]
    m.a[3,2] = a[5]
    m.a[4,0] = a[6]
    m.a[4,1] = a[7]
    m.a[4,2] = a[8]
    m.a[4,3] = a[9]
    m.a[5,0] = a[10]
    m.a[5,1] = a[11]
    m.a[5,2] = a[12]
    m.a[5,3] = a[13]
    m.a[5,4] = a[14]
    m.a[6,0] = a[15]
    m.a[6,1] = a[16]
    m.a[6,2] = a[17]
    m.a[6,3] = a[18]
    m.a[6,4] = a[19]
    m.a[6,5] = a[20]
    m.a[7,0] = a[21]
    m.a[7,1] = a[22]
    m.a[7,2] = a[23]
    m.a[7,3] = a[24]
    m.a[7,4] = a[25]
    m.a[7,5] = a[26]
    m.a[7,6] = a[27]
    m.a[8,0] = a[28]
    m.a[8,1] =  a[29]
    m.a[8,2]  = a[30]
    m.a[8,3] = a[31]
    m.a[8,4] = a[32]
    m.a[8,5] = a[33]
    m.a[8,6] = a[34]
    m.a[8,7] = a[35]
    m.a[9,0] =  a[36]
    m.a[9,1] = a[37]
    m.a[9,2] = a[38]
    m.a[9,3] = a[39]
    m.a[9,4] = a[40]
    m.a[9,5] =  a[41]
    m.a[9,6] = a[42]
    m.a[9,7] = a[43]
    m.a[9,8] = a[44]
    m.a[10,0] = a[45]
    m.a[10,1] = a[46]
    m.a[10,2] = a[47]
    m.a[10,3] = a[48]
    m.a[10,4] = a[49]
    m.a[10,5] = a[50]
    m.a[10,6] = a[51]
    m.a[10,7] = a[52]
    m.a[10,8] = a[53]
    m.a[10,9] = a[54]
    m.a[11,0] = a[55]
    m.a[11,1] = a[56]
    m.a[11,2] = a[57]
    m.a[11,3] = a[58]
    m.a[11,4] = a[59]
    m.a[11,5] = a[60]
    m.a[11,6] = a[61]
    m.a[11,7] = a[62]
    m.a[11,8] = a[63]
    m.a[11,9] = a[64]
    m.a[11,10]= a[65]
    m.a[12,0] = a[66]
    m.a[12,1] = a[67]
    m.a[12,2] = a[68]
    m.a[12,3] = a[69]
    m.a[12,4] = a[70]
    m.a[12,5] = a[71]
    m.a[12,6] = a[72]
    m.a[12,7] = a[73]
    m.a[12,8] = a[74]
    m.a[12,9] = a[75]
    m.a[12,10]= a[76]
    m.a[12,11]= a[77]
    m.a[13,0] = a[78]
    m.a[13,1] = a[79]
    m.a[13,2] = a[80]
    m.a[13,3] = a[81]
    m.a[13,4] = a[82]
    m.a[13,5] = a[83]
    m.a[13,6] = a[84]
    m.a[13,7] = a[85]
    m.a[13,8] = a[86]
    m.a[13,9] = a[87]
    m.a[13,10] =a[88]
    m.a[13,11] =a[89]
    m.a[13,12] =a[90]
    m.a[14,0] = a[91]
    m.a[14,1] = a[92]
    m.a[14,2] = a[93]
    m.a[14,3] = a[94]
    m.a[14,4] = a[95]
    m.a[14,5] = a[96]
    m.a[14,6] = a[97]
    m.a[14,7] = a[98]
    m.a[14,8] = a[99]
    m.a[14,9] = a[100]
    m.a[14,10] = a[101]
    m.a[14,11] = a[102]
    m.a[14,12] = a[103]
    m.a[14,13] = a[104]


 


    

     

    m.b[0] = b[0]
    m.b[1] = b[1]
    m.b[2] = b[2]
    m.b[3] = b[3]
    m.b[4] = b[4]
    m.b[5] = b[5]
    m.b[6] = b[6]
    m.b[7] = b[7]
    m.b[8] = b[8]
    m.b[9] = b[9]
    m.b[10] = b[10]
    m.b[11] = b[11]
    m.b[12] = b[12]
    m.b[13] = b[13]
    m.b[14] = b[14]

    y1 = +m.b[0]+m.b[1]+m.b[2]+m.b[3]+m.b[4]+m.b[5]+m.b[6]+m.b[7]+m.b[8]+m.b[9]+m.b[10]+m.b[11]+m.b[12]+m.b[13]+m.b[14] - (1/1)
    m.cons.add(y1 == 0 )

    y1 = +m.b[1]*c[1]+m.b[2]*c[2]+m.b[3]*c[3]+m.b[4]*c[4]+m.b[5]*c[5]+m.b[6]*c[6]+m.b[7]*c[7]+m.b[8]*c[8]+m.b[9]*c[9]+m.b[10]*c[10]+m.b[11]*c[11]+m.b[12]*c[12]+m.b[13]*c[13]+m.b[14]*c[14] - (1/2)
    m.cons.add(y1 == 0 )

    y1 = +m.b[1]*c[1]**2+m.b[2]*c[2]**2+m.b[3]*c[3]**2+m.b[4]*c[4]**2+m.b[5]*c[5]**2+m.b[6]*c[6]**2+m.b[7]*c[7]**2+m.b[8]*c[8]**2+m.b[9]*c[9]**2+m.b[10]*c[10]**2+m.b[11]*c[11]**2+m.b[12]*c[12]**2+m.b[13]*c[13]**2+m.b[14]*c[14]**2 - (1/3)
    m.cons.add(y1 == 0 )

    y1 = +m.b[2]*m.a[2,1]*c[1]+m.b[3]*m.a[3,1]*c[1]+m.b[3]*m.a[3,2]*c[2]+m.b[4]*m.a[4,1]*c[1]+m.b[4]*m.a[4,2]*c[2]+m.b[4]*m.a[4,3]*c[3]+m.b[5]*m.a[5,1]*c[1]+m.b[5]*m.a[5,2]*c[2]+m.b[5]*m.a[5,3]*c[3]+m.b[5]*m.a[5,4]*c[4]+m.b[6]*m.a[6,1]*c[1]+m.b[6]*m.a[6,2]*c[2]+m.b[6]*m.a[6,3]*c[3]+m.b[6]*m.a[6,4]*c[4]+m.b[6]*m.a[6,5]*c[5]+m.b[7]*m.a[7,1]*c[1]+m.b[7]*m.a[7,2]*c[2]+m.b[7]*m.a[7,3]*c[3]+m.b[7]*m.a[7,4]*c[4]+m.b[7]*m.a[7,5]*c[5]+m.b[7]*m.a[7,6]*c[6]+m.b[8]*m.a[8,1]*c[1]+m.b[8]*m.a[8,2]*c[2]+m.b[8]*m.a[8,3]*c[3]+m.b[8]*m.a[8,4]*c[4]+m.b[8]*m.a[8,5]*c[5]+m.b[8]*m.a[8,6]*c[6]+m.b[8]*m.a[8,7]*c[7]+m.b[9]*m.a[9,1]*c[1]+m.b[9]*m.a[9,2]*c[2]+m.b[9]*m.a[9,3]*c[3]+m.b[9]*m.a[9,4]*c[4]+m.b[9]*m.a[9,5]*c[5]+m.b[9]*m.a[9,6]*c[6]+m.b[9]*m.a[9,7]*c[7]+m.b[9]*m.a[9,8]*c[8]+m.b[10]*m.a[10,1]*c[1]+m.b[10]*m.a[10,2]*c[2]+m.b[10]*m.a[10,3]*c[3]+m.b[10]*m.a[10,4]*c[4]+m.b[10]*m.a[10,5]*c[5]+m.b[10]*m.a[10,6]*c[6]+m.b[10]*m.a[10,7]*c[7]+m.b[10]*m.a[10,8]*c[8]+m.b[10]*m.a[10,9]*c[9]+m.b[11]*m.a[11,1]*c[1]+m.b[11]*m.a[11,2]*c[2]+m.b[11]*m.a[11,3]*c[3]+m.b[11]*m.a[11,4]*c[4]+m.b[11]*m.a[11,5]*c[5]+m.b[11]*m.a[11,6]*c[6]+m.b[11]*m.a[11,7]*c[7]+m.b[11]*m.a[11,8]*c[8]+m.b[11]*m.a[11,9]*c[9]+m.b[11]*m.a[11,10]*c[10]+m.b[12]*m.a[12,1]*c[1]+m.b[12]*m.a[12,2]*c[2]+m.b[12]*m.a[12,3]*c[3]+m.b[12]*m.a[12,4]*c[4]+m.b[12]*m.a[12,5]*c[5]+m.b[12]*m.a[12,6]*c[6]+m.b[12]*m.a[12,7]*c[7]+m.b[12]*m.a[12,8]*c[8]+m.b[12]*m.a[12,9]*c[9]+m.b[12]*m.a[12,10]*c[10]+m.b[12]*m.a[12,11]*c[11]+m.b[13]*m.a[13,1]*c[1]+m.b[13]*m.a[13,2]*c[2]+m.b[13]*m.a[13,3]*c[3]+m.b[13]*m.a[13,4]*c[4]+m.b[13]*m.a[13,5]*c[5]+m.b[13]*m.a[13,6]*c[6]+m.b[13]*m.a[13,7]*c[7]+m.b[13]*m.a[13,8]*c[8]+m.b[13]*m.a[13,9]*c[9]+m.b[13]*m.a[13,10]*c[10]+m.b[13]*m.a[13,11]*c[11]+m.b[13]*m.a[13,12]*c[12]+m.b[14]*m.a[14,1]*c[1]+m.b[14]*m.a[14,2]*c[2]+m.b[14]*m.a[14,3]*c[3]+m.b[14]*m.a[14,4]*c[4]+m.b[14]*m.a[14,5]*c[5]+m.b[14]*m.a[14,6]*c[6]+m.b[14]*m.a[14,7]*c[7]+m.b[14]*m.a[14,8]*c[8]+m.b[14]*m.a[14,9]*c[9]+m.b[14]*m.a[14,10]*c[10]+m.b[14]*m.a[14,11]*c[11]+m.b[14]*m.a[14,12]*c[12]+m.b[14]*m.a[14,13]*c[13] - (1/6)
    m.cons.add(y1 == 0 )



#@title beta values { form-width: "100px" }

# used for linear order conditions plus linear stability polynomial matching conditions

import numpy as np

data=np.array([
    72.362092230726359621444062630453895436513878479204372239744118914339008361957274894021123680041692242,
1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
0.50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
0.166666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667,
0.0330278406028757527979332298466336725484247797530949148541780288765549097681738463834524909538164545,
0.00407886754615175430254075890828995614461800841500213845583064974804470591686586532729380258593128180,
0.000333527373172737242603705467098764792409315403431559178662909610705332831105918682590799425317677283,
0.00001884387445457719649926991124773774639255686232669238908076544455395077431328132833751461219355974457,
7.55565600673228187844472799772168548435456673497418681132813388498239377614064205083716425567987086e-7,
2.180674889659337806106474225333623253164152989661516940939984945939720123180192751320554488732457046e-8,
4.54329575207387209531051040170803945739047297526622495153899072480429169194765507413726877684905152e-10,
6.76964099198340395235741216662494203744813306696814739898497900746954960650250467206606180555417237e-12,
7.03301065055001143525874293209920770968902766712091555585888491039951221050037050929236960249104521e-14,
4.83760685817544104974918008895344303724730995701958486014397001982907347822621723661931383126126731e-16,
1.979536791110686622489750097380272593062776391994465218701741931323003217460574992708738821233104875e-18,
3.64746555383158031489459760645481664971731724579909139392166493612827917797499294769223945733650777e-21
              ])



#@title solve { form-width: "100px" }
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverFactory,SolverStatus,TerminationCondition
#from pyomo.core.expr import current as EXPR
import pyomo.core.expr as EXPR
from pyomo.core.expr import value
import time
from random import uniform,seed,randrange
import sys

reloadmodel=False

#data = np.loadtxt('beta.dat')
betaval = data[0]
betdic = dict(enumerate(data[1:]))

minbet=min(data)
#csos23 bettol=1e-9*minbet
#csos23
bettol=1e-6 ## Could push to 15?


m = pyo.ConcreteModel()


##NB EDIT THIS BY HAND
m.order = pyo.Param(initialize=3)
m.macc = pyo.Param(initialize=5)
print('Order:',m.order.value,'Acc:',m.macc.value)



#csos23 m.s = pyo.Param(initialize=len(betdic)-1)
m.s = pyo.Param(initialize=m.order.value*m.macc.value)


initval= 1/m.s.value
lneg=0 # Set to 0 to disable negative coefficients in bounding and initialisation: 0 appropriate for SSP
lim=1
skip=1

print('lim:',lim,'skip:',skip,'initval:',initval,'lneg:',lneg)

def setlsrk(m):
    for j in range(m.s.value-skip): #cols
        for i in range(j+skip+1,m.s.value):
            m.cons.add(m.a[i,j] == m.b[j] )


def azeros(m):
    for i in m.rows:
        for j in range(i,m.s.value): #cols
            m.cons.add(m.a[i,j] == 0 )


m.rows= range(m.s.value)
m.cols= range(m.s.value)

m.a= pyo.Var( m.rows, m.cols, bounds=(-lneg*lim,lim))
m.b= pyo.Var( m.rows, bounds=(-lneg*lim,lim))
m.beta=pyo.Param(initialize=betaval)  #beta_R
m.bet = pyo.Param(range(m.s.value+1), initialize=betdic) #beta_j

print('beta_R:',m.beta.value)

m.cons = pyo.ConstraintList()

mycons(m)
#csos23
setlsrk(m)
#csos23?
#azeros(m) #Nones faster
mybets(m)

#m.pprint()


seedval = randrange(sys.maxsize)
#seedval=2842383511757470787
seed(seedval)
print('seedval:',seedval)



### solver = SolverFactory('ipopt',executable='/content/Ipopt-3.12.13/build/bin/ipopt')
solver = SolverFactory('ipopt',executable='ipopt')
#solver.options['max_iter'] = 10000
solver.options['nlp_scaling_method'] = 'none'
solver.options['linear_solver'] = 'ma27'
solver.options['honor_original_bounds']='yes'

solver.options['bound_relax_factor']=bettol ## default 1e-8   https://kiwi.oden.utexas.edu/papers/Multi-output-multi-fidelity-MLBLUE-Croci-Wright-Willcox.pdf

solver.options['tol'] = bettol # default 1e-8
solver.options['dual_inf_tol']=1e8*bettol #default 1, Larger more accurate
solver.options['constr_viol_tol']=1e4*bettol #default 1e-4
solver.options['compl_inf_tol']=1e4*bettol #default 1e-4

solver.options['acceptable_tol'] = bettol # default 1e-6
solver.options['acceptable_dual_inf_tol']=1e8*bettol #default 1e10
solver.options['acceptable_constr_viol_tol']=1e4*bettol #default 1e-2
solver.options['acceptable_compl_inf_tol']=1e4*bettol #default 1e-2



niter=1

#solver = SolverFactory('multistart')
#niter=1
#nitermulti=10

#solver = SolverFactory('couenne')
#niter=1


if reloadmodel:
    with open('m.pkl', mode='rb') as file: #reload model
        m = cloudpickle.load(file)
        print('Reload')
else:
    objvalbest=1e10
    for iter in range(niter):

        start = time.time()

        for i in range(1,m.s.value):
            for j in range(i):
                #csos23
                m.a[i,j].value=uniform(-lneg*initval,initval)
                #m.a[i,j].value=initval
        for j in m.cols:
            #csos23
            m.b[i].value=uniform(-lneg*initval,initval)
            #m.b[i].value=initval   ## sum b=1 !



        try:

            results=solver.solve(m,tee=True)
            #results=solver.solve(m,tee=False)

            #results=solver.solve(m,solver='ipopt',strategy='rand_distributed',iterations=nitermulti)
            if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
                end = time.time()
                objval=m.obj()
                print(iter,'Time:',end - start,'Objective:',objval)
                if objval < objvalbest:
                    objvalbest=objval
                    mbest=m.clone()
                    resultsbest=results
            elif (results.solver.termination_condition == TerminationCondition.infeasible):
                print(iter,'Infeasible.')
            else:
                print(iter,'Failed. Solver Status:',  result.solver.status)
        except:
            print(iter,'Ipopt error, moving on.')

        #m.display()
        #results.write()


    m=mbest.clone()  # revert to best instance
    results=resultsbest


    import cloudpickle  # write model to file
    with open('m.pkl', mode='wb') as file:
        cloudpickle.dump(m, file)


objval=m.obj()
print('Final objective value:',objval)

#@title stability polynomials { form-width: "100px" }


a = []
for i in range(1,m.s.value):
    for j in range(i):
              
            c = ('%5.20f'%m.a[i,j].value)
            a.append(float(c))

#print("\n")
b = []
#print("B")
for j in m.cols:
    c = ('%5.20f'%m.b[j].value)
    b.append(float(c))

print(a)
print(b)

	
