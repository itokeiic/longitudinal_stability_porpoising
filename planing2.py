from lanczos_approximation import *
from scipy.optimize import *
from scipy.linalg import *
from scipy import *
import pylab
pylab.rc('text',usetex=True)
pylab.rcParams['figure.subplot.bottom']=0.12
pylab.rcParams['figure.subplot.left']=0.16
#from math import *
def trim(M,U,B,beta_deg,l_cg,vcg,f,epsilon,rho=1026.0,g=9.81,AHR=150.0E-6,nu=1.19E-6):
    a=vcg-B/4.0*tan(beta_deg*pi/180)
    maxiter=100000
    def trim_angle(tau_deg,M):
        F_n_B=U/sqrt(g*B)
        tau=tau_deg*pi/180.0
        #tau_deg=tau*180/pi
        C_L_beta=M*g/(.5*rho*U**2*B**2)
        #C_L_0=linspace(0,1);
        def f1(C_L_0):
            y=C_L_0-.0065*beta_deg*C_L_0**.60-C_L_beta
            return y
        C_L_0=brentq(f1,0.0,2.0,maxiter=maxiter)
        
        def f2(lambda_w):
            y=tau_deg**1.1*(.012*lambda_w**.5+.0055*lambda_w**2.5/F_n_B**2)-C_L_0
            return y
        lambda_w=brentq(f2,0.0,500.0,maxiter=maxiter)
        if lambda_w > 4.0:
            print 'lambda_w is grater than 4, and is out of the range of validity.' 
        x_s=B*tan(beta_deg*pi/180)/(pi*tau)
        L_C=B*lambda_w-.5*x_s
        if L_C < 0.0:
            print 'chine not wet!'
            L_C=0.0
        L_K=2*B*lambda_w-L_C
        D=L_K*sin(tau)
        
        S1=pi/2*tau*x_s**2/sin(beta_deg*pi/180)
        S2=B/cos(beta_deg*pi/180)*L_C
        S=S1+S2
        R_n=U*L_K/nu
        C_F_ITTC=.075/(log10(R_n)-2)**2
        DeltaC_F=(44*((AHR/L_K)**(1.0/3.0)-10*R_n**(-1.0/3.0))+.125)*10**(-3)
        C_F=C_F_ITTC+DeltaC_F
        R_v=.5*rho*C_F*S*U**2
        R_p=M*g*tan(tau)
        Drag_hydro=R_v*cos(tau)+R_p
        
        
        l_p=(.75-1/(5.21*(F_n_B**2/lambda_w**2)+2.39))*lambda_w*B
        c=l_cg-l_p;
        res=c/cos(tau)*(M*g-(M*g*sin(tau)+R_v)/cos(epsilon)*sin(tau+epsilon)+R_v*sin(tau))+R_v*a-(M*g*sin(tau)+R_v)/cos(epsilon)*f
        return res,F_n_B,C_L_beta,C_L_0,lambda_w,x_s,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p
    
    def func(tau_deg,M):
        res,F_n_B,C_L_beta,C_L_0,lambda_w,x_s,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=trim_angle(tau_deg,M)
        return res
    # lowtrim=0.2;hightrim=45.5
    lowtrim=0.1;hightrim=45.5
    print 'trim residue low, high:',func(lowtrim,M),func(hightrim,M)
    
    tau_deg=brentq(func,lowtrim,hightrim,xtol=1.0e-6,maxiter=maxiter,args=(M))
    res,F_n_B,C_L_beta,C_L_0,lambda_w,x_s,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=trim_angle(tau_deg,M)
    if L_C < 0.0:
        print 'chine not wet. setting L_C = 0'
        L_C=0.0
    else:
        print 'wet chine!'
    if lambda_w > 4.0:
        print 'lambda_w is grater than 4, and is out of the range of validity.' 
    else:
        print 'lambda_w is within range of validity!'

    return (res,tau_deg,F_n_B,C_L_beta,C_L_0,lambda_w,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p)


#test values======================================
# rho=1026.0
# g=9.81
# beta_deg=20.0
# tau_deg=4.0
# beta=beta_deg*pi/180.0
# tau=tau_deg*pi/180.0
# lambda_w=4.0
# B=4.0
# vcg=0.25*B
# l_cg=2.13*B
# M=1.28*rho*B**3
# I55=M*(1.3*B)**2
# F_n_B=5.0
# l_p=(.75-1.0/(5.21*(F_n_B/lambda_w)**2+2.39))*lambda_w*B
# #l_p=2.0/3.0*lambda_w*B
# x_s=B*tan(beta_deg*pi/180)/(pi*tau)
# L_C=B*lambda_w-.5*x_s
# L_K=2*B*lambda_w-L_C
# U=F_n_B*sqrt(g*B)
# print 'U = ',U
# C_L_beta=M*g/(.5*rho*U**2*B**2)
# print 'C_L_beta = ',C_L_beta
# def f1(C_L_0):
    # y=C_L_0-.0065*beta_deg*C_L_0**.60-C_L_beta
    # return y
# C_L_0=brentq(f1,0.0,5.0)
# print 'C_L_0 = ',C_L_0
# C_L_0=tau_deg**1.1*(0.012*lambda_w**0.5+0.0055*lambda_w**2.5/F_n_B**2)
# print 'C_L_0 = ',C_L_0
# C_L_beta=C_L_0-0.0065*beta_deg*C_L_0**0.60 #2015-05-12 this makes c53,c55 values close to Faltinsen's book's value
# print 'C_L_beta = ',C_L_beta
# D=L_K*sin(tau)
#================================================
def porpoising(U,M,I55,B,beta_deg,tau_deg,vcg,L_K,L_C,l_cg,l_p,C_L_0,C_L_beta,lambda_w,rho=1026.0,g=9.81):
    F_n_B=U/sqrt(g*B)
    x_s=L_K-L_C
    beta=beta_deg*pi/180.0
    tau=tau_deg*pi/180.0
    b=[0, 4, 7.5, 10, 15, 20, 25, 30, 40]
    z=[pi/2-1, .5695, .5623, .5556, .5361, .5087, .4709, .4243, .2866]
    # one_p_zmax_o_Vt=pi/2.0
    one_p_zmax_o_Vt=1.0+interp(beta_deg,b,z)
    zwl=vcg*cos(tau)-(L_K-l_cg)*sin(tau)
    print 'zwl = ',zwl
    dlambda_w_o_deta_5=(-vcg/B)/(sin(tau))**2+(zwl/B)/(sin(tau))**2*cos(tau)+0.25*tan(beta)/(one_p_zmax_o_Vt*tau*tau)
    dlambda_w_o_deta_3=-1.0/(sin(tau)*B)
    print 'dlambda_w/deta_5 = ',dlambda_w_o_deta_5
    print 'dlambda_w/deta_3 = ',dlambda_w_o_deta_3
    dC_L_0_o_deta_5=1.1*(180.0/pi)**1.1*tau**0.1*(0.012*lambda_w**0.5+0.0055*lambda_w**2.5/F_n_B**2)+tau_deg**1.1*(0.006*lambda_w**(-0.5)+0.01375*lambda_w**1.5/F_n_B**2)*dlambda_w_o_deta_5
    dC_L_0_o_deta_3=tau_deg**1.1*(0.006*lambda_w**(-0.5)+0.01375*lambda_w**1.5/F_n_B**2)*dlambda_w_o_deta_3
    print 'dC_L_0/deta_5 = ',dC_L_0_o_deta_5
    print 'dC_L_0/deta_3 = ',dC_L_0_o_deta_3
    print 'B:',B,'beta_deg:',beta_deg,'C_L_0:',C_L_0
    c33=-B*dC_L_0_o_deta_3*(1-0.0039*beta_deg*C_L_0**(-0.4))
    C33=c33*(0.5*rho*U*U*B)
    c35=-dC_L_0_o_deta_5*(1-0.0039*beta_deg*C_L_0**(-0.4))
    C35=c35*(0.5*rho*U*U*B*B)
    print 'c33 = ',c33*0.5*rho*U*U*B/(rho*g*B*B)
    print 'c35 = ',c35*0.5*rho*U*U*B*B/(rho*g*B*B*B)
    dl_p_o_dlambda_w=(0.75-(15.63*(F_n_B**2/lambda_w**2)+2.39)/(5.21*(F_n_B**2/lambda_w**2)+2.39)**2)*B
    c53=-(dl_p_o_dlambda_w*dlambda_w_o_deta_3*C_L_beta+(l_p-l_cg)/B*(-c33))
    C53=c53*(0.5*rho*U*U*B*B)
    c55=-(dl_p_o_dlambda_w/B*dlambda_w_o_deta_5*C_L_beta+(l_p-l_cg)/B*(-c35))
    C55=c55*(0.5*rho*U*U*B*B*B)
    print 'c53 = ',c53*0.5*rho*U*U*B*B/(rho*g*B*B*B)
    print 'c55 = ',c55*0.5*rho*U*U*B*B*B/(rho*g*B*B*B*B)
    
    K=(1.0/tan(beta))*((pi/sin(beta))*gamma(1.5-beta/pi)/(gamma(1.0-beta/pi)**2*gamma(0.5+beta/pi))-1.0)
    print 'K = ',K.real
    a33_1=(K.real/24.0)*tan(beta)**3/(one_p_zmax_o_Vt*tau)
    print 'a33_1 = ',a33_1
    x_G=L_K-l_cg
    print 'x_G = ',x_G
    a35_1=a53_1=a33_1*x_G/B-(K.real/64.0)*tan(beta)**4/(one_p_zmax_o_Vt*tau)**2
    print 'a35_1 = ',a35_1,'a53_1 = ',a53_1
    a55_1=(K.real/160.0)*tan(beta)**5/(one_p_zmax_o_Vt*tau)**3-(K.real/32.0)*(x_G/B)*tan(beta)**4/(one_p_zmax_o_Vt*tau)**2+(x_G/B)**2*a33_1
    print 'a55_1 = ',a55_1
    C_1=(2*tan(beta)**2/pi)*K.real
    a33_2=C_1*(pi/8.0)*L_C/B
    print 'a33_2 = ',a33_2
    a35_2=a53_2=-C_1*pi/16.0*((L_K/B)**2-(x_s/B)**2)+(x_G/B)*a33_2
    print 'a35_2 = ',a35_2,'a53_2 = ',a53_2
    a55_2=C_1*pi/24.0*((L_K/B)**3-(x_s/B)**3)-C_1*pi/8.0*(x_G/B)*((L_K/B)**2-(x_s/B)**2)+(x_G/B)**2*a33_2
    print 'a55_2 = ',a55_2
    print 'a33 = ',a33_1+a33_2
    print 'a35 = ',a35_1+a35_2
    print 'a53 = ',a53_1+a53_2
    print 'a55 = ',a55_1+a55_2
    A33=(a33_1+a33_2)*(rho*B**3)
    A35=(a35_1+a35_2)*(rho*B**4)
    A53=(a53_1+a53_2)*(rho*B**4)
    A55=(a55_1+a55_2)*(rho*B**5)
    C_L_0=tau_deg**1.1*0.012*lambda_w**0.5 #modified C_L_0 assuming F_n_B=inf
    print 'C_L_0 = ',C_L_0
    dC_L_0_o_dtau=(180.0/pi)**1.1*0.0132*tau**0.1*lambda_w**0.5
    dC_L_beta_o_dtau=dC_L_0_o_dtau*(1.0-0.0039*beta_deg*C_L_0**(-0.4))
    b33=0.5*F_n_B*dC_L_beta_o_dtau
    B33=b33*(rho*B**3*sqrt(g/B))
    print 'b33 = ',b33
    b53=0.75*lambda_w-l_cg/B
    print 'b53 = ',b53
    B53=b53*(B33*B)
    print 'b53 = ',B53/(rho*B**4*sqrt(g/B))
    d=0.5*B*tan(beta)
    #d=L_K*sin(tau)
    a33sec=rho*d*d*K.real
    #x_T=l_cg*cos(tau)-vcg*sin(tau) 
    x_T=l_cg #2015-05-12 this makes b55 values close to Faltinsen's book's value
    B55=U*x_T**2*a33sec
    print 'b55 = ',B55/(rho*B**5*sqrt(g/B))
    B35=-U*A33-U*x_T*a33sec
    print 'b35 = ',B35/(rho*B**4*sqrt(g/B))
    AA=zeros((4,4))
    BB=zeros((4,4))
    AA[0,0]=M+A33
    AA[0,1]=A35
    AA[1,0]=A53
    AA[1,1]=I55+A55
    AA[2,2]=1.0
    AA[3,3]=1.0
    BB[0,0]=B33
    BB[0,1]=B35
    BB[0,2]=C33
    BB[0,3]=C35
    BB[1,0]=B53
    BB[1,1]=B55
    BB[1,2]=C53
    BB[1,3]=C55
    BB[2,0]=-1.0
    BB[3,1]=-1.0
    print 'AA = ',AA
    print 'BB = ',BB
    KK=dot(inv(AA),-BB)
    print 'KK = ',KK
    w,vr=eig(KK)
    # print 'eigenvalues:',w
    # print 'maximum eigen value = ',max(w)
    # print 'real part of maximum eigen value = ',max(w).real 
    return (w,vr,KK,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55)
    
def main():
    figformat='pdf'
    figdpi=600
    #input data for towing tank model===================================
    plotresol=20
    rho=1000.0 #water density
    g=9.81
    M=8.09#27000 #weight of seaplane in kg
    U=5#20.58+0.0 #planing speed in m/s
    B=0.2 # 4.27 #beam length in m
    I55=M*0.608**2#M*(1.3*B)**2
    beta_deg=8.0 #deadrise angle in degrees
    l_cg=1.54/2.0-0.666#8.84+0.0 #longitudinal location of center of gravity from step
    vcg=0.453#B/4.0*tan(beta_deg*pi/180)+8.
    # a=vcg-B/4.0*tan(beta_deg*pi/180)
    f=0.0
    epsilon=0.0*pi/180 #thrust angle with respect to keel line
    AHR=150E-6 #average hull roughness
    nu=1.52E-6 #kinetic viscosity of water
    #===================================================================
    #input data for Faltinsen's trim example=================================
    # plotresol=2
    # rho=1026.0 #saltwater density
    # g=9.81
    # M=27000 #weight of seaplane in kg
    # U=20.58+0.0 #planing speed in m/s
    # B=4.27 #beam length in m
    # I55=M*(1.3*B)**2
    # beta_deg=10.0 #deadrise angle in degrees
    # l_cg=8.84+0.0 #longitudinal location of center of gravity from step
    # vcg=B/4.0*tan(beta_deg*pi/180)#+0.42
    # a=vcg-B/4.0*tan(beta_deg*pi/180)
    # f=0#0.15
    # epsilon=0#4.0*pi/180 #thrust angle with respect to keel line
    # AHR=150E-6 #average hull roughness
    # nu=1.19E-6 #kinetic viscosity of saltwater
    #===================================================================
    #input data for Faltinsen's porposing example=================================
    # plotresol=2
    # rho=1026.0 #saltwater density
    # g=9.81
    # tau_deg=4.0
    # B=4.0 #beam length in m
    # M=1.28*rho*B**3 #weight of seaplane in kg
    # U=5.0*(g*B)**0.5 #planing speed in m/s
    # I55=M*(1.3*B)**2
    # beta_deg=20.0 #deadrise angle in degrees
    # l_cg=2.13*B #longitudinal location of center of gravity from step
    # vcg=0.25*B #B/4.0*tan(beta_deg*pi/180)#+0.42
    # a=vcg-B/4.0*tan(beta_deg*pi/180)
    # lambda_w = 4.0
    # f=0#0.15
    # epsilon=0#4.0*pi/180 #thrust angle with respect to keel line
    # AHR=150E-6 #average hull roughness
    # nu=1.19E-6 #kinetic viscosity of saltwater
    #===================================================================

    res,tau_deg,F_n_B,C_L_beta,C_L_0,lambda_w,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=trim(M,U,B,beta_deg,l_cg,vcg,f,epsilon,rho=rho,g=g,AHR=AHR,nu=nu)
    
    print "Trim angle, tau = ",tau_deg,"(deg)\n", 
    print "Froude Number with respect to the beam, F_n_B = ",F_n_B
    print "Lift Coefficient at given deadrise angle, C_L_beta = ", C_L_beta
    print "Equivalent lift coefficient with zero deadrise angle, C_L_0 = ",C_L_0
    print "Mean wetted length to beam ratio, lambda_w = ",lambda_w
    print "The difference in wetted length L_K-L_C, x_s = ", L_K-L_C
    print "Chine wetted length, L_C = ", L_C
    print "Keel wetted length, L_K = ", L_K
    print "Draft, D = ", D
    print "Wetted area, S = ", S
    print "Reynolds Number with respect to L_K, R_n = ", R_n
    print "Friction coefficient, C_F = ", C_F
    print "Longitudinal frictional force, R_v = ", R_v
    print "Horizontal component of hydrodynamic pressure, R_p = ", R_p
    print "Hydrodynamic Drag = ", Drag_hydro
    print "Required Power = %f kW or %f hp"%(Drag_hydro*U/1000,(Drag_hydro*U/1000)*1.359622)
    print "Centre of pressure measured from the transom along the keel, l_p = ",l_p
    print "residue = ", res
    print ' '
    w,vr,KK,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55= porpoising(U,M,I55,B,beta_deg,tau_deg,vcg,L_K,L_C,l_cg,l_p,C_L_0,C_L_beta,lambda_w,rho=rho,g=g)
    print ' '
    print 'eigenvalues:\n',w
    print 'eigenvectors:\n',vr
    print 'real part of maximum eigen value = ',max(w).real  
    
    fignum=0
    fignum=fignum+1
    pylab.figure(fignum)
    recxset=[];recmaxwset=[];reclcg=[]
    # for l_cg,I55 in zip([1.54/2.0-0.666,1.54/2.0-0.764],[M*0.608**2,M*0.799**2]):
    for l_cg in [1.54/2.0-0.64,1.54/2.0-0.67]:
        reclcg.append(l_cg)
        recx=[];recmaxw=[]
        for U in linspace(4.0,8.0,plotresol):
            
            res,tau_deg,F_n_B,C_L_beta,C_L_0,lambda_w,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=trim(M,U,B,beta_deg,l_cg,vcg,f,epsilon,rho=rho,g=g,AHR=AHR,nu=nu)
            w,vr,KK,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55= porpoising(U,M,I55,B,beta_deg,tau_deg,vcg,L_K,L_C,l_cg,l_p,C_L_0,C_L_beta,lambda_w,rho=rho,g=g)
            recx.append(F_n_B); recmaxw.append(max(w).real*sqrt(B/g))
            print '**moment residue = ',res
            print '**max_eig_real*sqrt(B/g) = ',max(w).real*sqrt(B/g)
        recxset.append(list(recx));recmaxwset.append(list(recmaxw))
    pylab.plot(recxset[0],recmaxwset[0],recxset[1],recmaxwset[1],'--')
    pylab.legend((r'$l_{cg}/B=$'+str(reclcg[0]/B),r'$l_{cg}/B=$'+str(reclcg[1]/B)),loc=0)
    pylab.grid()
    pylab.xlabel(r'$F_{n_B}$',fontsize=18)
    pylab.ylabel('$\Re(\sigma)_{\max}\sqrt{B/g}$',fontsize=18)   
    pylab.xticks(fontsize=18)
    pylab.yticks(fontsize=18)
    pylab.legend(loc=0,prop={'size':12})
    
    pylab.savefig('eig_vs_Fn.%s'%figformat,format=figformat,dpi=figdpi)
    l_cg=1.54/2.0-0.666
    fignum=fignum+1
    pylab.figure(fignum)
    # for i,(l_cg,I55) in enumerate(zip([1.54/2.0-0.666,1.54/2.0-0.764],[M*0.608**2,M*0.799**2])):
    for i,l_cg in enumerate([1.54/2.0-0.64,1.54/2.0-0.67]):
        recx=[];recmaxw=[]
        for U in linspace(4.0,8.0,plotresol):
            res,tau_deg,F_n_B,C_L_beta,C_L_0,lambda_w,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=trim(M,U,B,beta_deg,l_cg,vcg,f,epsilon,rho=rho,g=g,AHR=AHR,nu=nu)
            recx.append(F_n_B); recmaxw.append(tau_deg)
            print '**moment residue = ',res
            print '**max_eig_real*sqrt(B/g) = ',max(w).real*sqrt(B/g)

        if i==0:
            pylab.plot(recx,recmaxw,label=r'$l_{cg}/B=$'+str(l_cg/B))
        if i==1:
            pylab.plot(recx,recmaxw,'--',label=r'$l_{cg}/B=$'+str(l_cg/B))
    pylab.legend(loc=0)
    pylab.grid()
    pylab.xlabel(r'$F_{n_B}$',fontsize=18)
    pylab.ylabel(r'$\tau$ [deg]',fontsize=18)    
    pylab.legend()
    pylab.xticks(fontsize=18)
    pylab.yticks(fontsize=18)
    pylab.legend(loc=0,prop={'size':12})

    pylab.savefig('tau_vs_Fn.%s'%figformat,format=figformat,dpi=figdpi)
    # U=6.0
    xlcg=linspace(0.03,1.0,plotresol)
    yvcg=linspace(0.001,1.0,plotresol)
    X,Y=pylab.meshgrid(xlcg,yvcg)
    for U in [4.0, 6.0, 8.0]:
        fignum=fignum+1
        pylab.figure(fignum)
        recmaxw=[]
        for l_cg,vcg in zip(X.flatten(),Y.flatten()):
            print "U = %f, l_cg = %f, vcg = %f"%(U,l_cg,vcg)
            res,tau_deg,F_n_B,C_L_beta,C_L_0,lambda_w,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=trim(M,U,B,beta_deg,l_cg,vcg,f,epsilon,rho=1000.0,g=g,AHR=AHR,nu=nu)
            w,vr,KK,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55= porpoising(U,M,I55,B,beta_deg,tau_deg,vcg,L_K,L_C,l_cg,l_p,C_L_0,C_L_beta,lambda_w,rho=rho,g=g)
            recmaxw.append(max(w).real*sqrt(B/g))
            print '**moment residue = ',res
            print '**max_eig_real*sqrt(B/g) = ',max(w).real*sqrt(B/g)
        recmaxw=array(recmaxw)
        Z=recmaxw.reshape(len(xlcg),len(yvcg))
        cs=pylab.contour(X/B,Y/B,Z,10)
        pylab.clabel(cs, inline=1, fontsize=10)
        pylab.grid()
        pylab.xlabel(r'$l_{cg}/B$',fontsize=18)
        pylab.ylabel(r'$vcg/B$',fontsize=18)
        pylab.title(r'Maximum Real Part of Eigenvalues, $\Re(\sigma)_{\max}\sqrt{B/g}$')
        pylab.xticks(fontsize=18)
        pylab.yticks(fontsize=18)

        pylab.savefig('contour_vcg_vs_lcg_U%d.%s'%(U,figformat),format=figformat,dpi=figdpi)
    pylab.show()

if __name__ == "__main__":
    main()
