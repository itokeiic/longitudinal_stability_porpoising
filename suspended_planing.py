from math import sqrt
import scipy
import scipy.linalg
import pylab
import planing2
pylab.rc('text',usetex=True)
pylab.rcParams['figure.subplot.bottom']=0.12
pylab.rcParams['figure.subplot.left']=0.16
def Maxi(a):
    for l,v in enumerate(a):
       if l==0:
          r=v
          index=l
       else:
          if r<v:
             r=v
             index=l
    return [r, index]


def susporpoising(m_A,I_A,m_B,I_B,kf,kb,cf,cb,l_Af,l_Ab,l_Bf,l_Bb,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55):
    A=scipy.zeros((8,8))
    B=scipy.zeros((8,8))

    A[0,0]=m_A
    A[1,1]=I_A
    A[2,2]=1.
    A[3,3]=1.
    A[4,4]=m_B+A33
    A[4,5]=A35
    A[5,4]=A53
    A[5,5]=I_B+A55
    A[6,6]=1.
    A[7,7]=1.
    
    B[0,0]=cf+cb
    B[0,1]=cf*l_Af-cb*l_Ab
    B[0,2]=kf+kb
    B[0,3]=kf*l_Af-kb*l_Ab
    B[0,4]=-cf-cb
    B[0,5]=-cf*l_Bf+cb*l_Bb
    B[0,6]=-kf-kb
    B[0,7]=-kf*l_Bf+kb*l_Bb
    B[1,0]=cf*l_Af-cb*l_Ab
    B[1,1]=cf*l_Af**2+cb*l_Ab**2
    B[1,2]=kf*l_Af-kb*l_Ab
    B[1,3]=kf*l_Af**2+kb*l_Ab**2
    B[1,4]=-cf*l_Af+cb*l_Ab
    B[1,5]=-cf*l_Af*l_Bf-cb*l_Ab*l_Bb
    B[1,6]=-kf*l_Af+kb*l_Ab
    B[1,7]=-kf*l_Af*l_Bf-kb*l_Ab*l_Bb
    B[4,0]=-cf-cb
    B[4,1]=-cf*l_Af+cb*l_Ab
    B[4,2]=-kf-kb
    B[4,3]=-kf*l_Af+kb*l_Ab
    B[4,4]=B33+cf+cb
    B[4,5]=B35+cf*l_Bf-cb*l_Bb
    B[4,6]=C33+kf+kb
    B[4,7]=C35+kf*l_Bf-kb*l_Bb
    B[5,0]=-cf*l_Bf+cb*l_Bb
    B[5,1]=-cf*l_Af*l_Bf-cb*l_Ab*l_Bb
    B[5,2]=-kf*l_Bf+kb*l_Bb
    B[5,3]=-kf*l_Af*l_Bf-kb*l_Ab*l_Bb
    B[5,4]=B53+cf*l_Bf-cb*l_Bb
    B[5,5]=B55+cf*l_Bf**2+cb*l_Bb**2
    B[5,6]=C53+kf*l_Bf-kb*l_Bb
    B[5,7]=C55+kf*l_Bf**2+kb*l_Bb**2
    B[2,0]=-1.
    B[3,1]=-1.
    B[6,4]=-1.
    B[7,5]=-1.
    K=scipy.dot(scipy.linalg.inv(A),-B)
    ws,vrs=scipy.linalg.eig(K)
    return (ws,vrs,K)
    
def main():
    #input data=========================================================
    rho=1000.0 #saltwater density
    g=9.81
    M=8.09#27000 #weight of seaplane in kg
    U=7.0#20.58+0.0 #planing speed in m/s
    B=0.2 # 4.27 #beam length in m
    I55=M*0.608**2#M*(1.3*B)**2
    beta_deg=8.0 #deadrise angle in degrees
    l_cg=1.54/2.0-0.66#8.84+0.0 #longitudinal location of center of gravity from step
    vcg=0.453#B/4.0*tan(beta_deg*pi/180)+8.
    a=vcg-B/4.0*scipy.tan(beta_deg*scipy.pi/180)
    f=0.0
    epsilon=0.0*scipy.pi/180 #thrust angle with respect to keel line
    AHR=150E-6 #average hull roughness
    nu=1.19E-6 #kinetic viscosity of saltwater
    #===================================================================
    res,tau_deg,F_n_B,C_L_beta,C_L_0,lambda_w,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=planing2.trim(M,U,B,beta_deg,l_cg,vcg,f,epsilon,rho=1000.0,g=9.81,AHR=150.0E-6,nu=1.19E-6)
    
    w,vr,KK,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55= planing2.porpoising(U,M,I55,B,beta_deg,tau_deg,vcg,L_K,L_C,l_cg,l_p,C_L_0,C_L_beta,lambda_w,rho=1000.0,g=9.81)
    
    m_A=10.79/2
    #I_A=m_A*0.61**2+0.3357
    m_B=5.39/2
    I_B=m_B*0.55**2
    I_A=I55-I_B
    cf=10.
    cb=10.
    kf=2117.
    kb=2117.
    l_Af=0.2
    l_Ab=0.2
    l_Bf=0.2
    l_Bb=0.2
    ws,vrs,K=susporpoising(m_A,I_A,m_B,I_B,kf,kb,cf,cb,l_Af,l_Ab,l_Bf,l_Bb,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55)
    print 'K=\n',K
    #w,vr=scipy.linalg.eig(KK)
    
    print 'eigenvalues of rigid case: \n',w
    print 'real part of maximum eigen value = ',max(w).real 
     
    print 'eigenvalues:\n',ws
    print 'real part of maximum eigen value = ',max(ws).real
    
    kf_list=[1.0e12,2117.0,2117.0,2117.0,2117.0,2117.0,200.0,1.0e12,2117.0,2117.0,2117.0]
    kb_list=[1.0e12,2117.0,2117.0,2117.0,2117.0,2117.0,200.0,2117.0,1.0e12,2117.0,2117.0]
    cf_list=[10.0,0.0,20.0,30.0,50.0,80.0,0.0,0.0,0.0,80.0,0.0]
    cb_list=[10.0,0.0,20.0,30.0,50.0,80.0,0.0,0.0,0.0,0.0,80.0]
    # kf_list=[1.0e12]
    # kb_list=[1.0e12]
    # cf_list=[10]
    # cb_list=[10]

    for i,(kf,kb,cf,cb) in enumerate(zip(kf_list,kb_list,cf_list,cb_list)):
        x=[];y=[];susy=[]
        for U in scipy.linspace(4.0,8.5,20):
            res,tau_deg,F_n_B,C_L_beta,C_L_0,lambda_w,L_C,L_K,D,S,R_n,C_F, R_v,R_p,Drag_hydro,l_p=planing2.trim(M,U,B,beta_deg,l_cg,vcg,f,epsilon,rho=1000.0,g=9.81,AHR=150.0E-6,nu=1.19E-6)
            
            w,vr,KK,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55= planing2.porpoising(U,M,I55,B,beta_deg,tau_deg,vcg,L_K,L_C,l_cg,l_p,C_L_0,C_L_beta,lambda_w,rho=1000.0,g=9.81)
            print 'rigid eigenvalues:\n',w
            print 'rigid eigenvectors:\n',vr
            ws,vrs,K=susporpoising(m_A,I_A,m_B,I_B,kf,kb,cf,cb,l_Af,l_Ab,l_Bf,l_Bb,A33,A35,A53,A55,B33,B35,B53,B55,C33,C35,C53,C55)
            print 'flexible eigenvalues:\n',ws
            print 'flexible eigenvectors:\n',vrs
            x.append(F_n_B);y.append(max(w).real*sqrt(B/g));susy.append(max(ws).real*sqrt(B/g))
        pylab.figure(i+1)
        pylab.plot(x,y,label=r'Rigid support')
        pylab.plot(x,susy,'--',label=r'$k_f B/(Mg)=$'+'%5.3e'%(kf*B/(M*g))+r', $k_b B/(Mg)=$'+'%5.3e'%(kb*B/(M*g))+',\n'+ r'$c_f \sqrt{B/g}/M=$'+'%5.3e'%(cf*sqrt(B/g)/M)+r', $c_b \sqrt{B/g}/M=$'+'%5.3e'%(cb*sqrt(B/g)/M))
        pylab.xlabel(r'$F_{n_B}$',fontsize=18)
        pylab.ylabel(r'$\Re(\sigma)_{\max}\sqrt{B/g}$',fontsize=18)
        pylab.xticks(fontsize=18)
        pylab.yticks(fontsize=18)
        pylab.legend(loc=0,prop={'size':12})
        pylab.grid()
        pylab.savefig('flexible_support_%d.pdf'%i,format='pdf')
    pylab.show()
if __name__ == "__main__":
    main()
