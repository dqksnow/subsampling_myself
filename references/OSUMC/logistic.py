#!/usr/bin/env python
# coding: utf-8

import numpy as np
from sklearn.linear_model.logistic import LogisticRegression
from scipy.optimize import fsolve
from scipy.optimize import fmin_ncg
from scipy.optimize import fmin_cg
from scipy.linalg import solve
import scipy.sparse as sp
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import numba as nb
import time
get_ipython().run_line_magic('matplotlib', 'inline')
warnings.filterwarnings("ignore")


N = 100000
#N=25000
d = 20
Sigma = 0.5 * np.ones((d,d), dtype=float)
for i in range(d):
    Sigma[i,i] = 1
mu_mz = np.zeros(d)
mu_nz = 0.5 * np.ones(d)
U1 = np.diag(1/np.array(range(1,d+1)))
Sigma1 = np.dot(np.dot(U1, Sigma), U1)
beta0 = np.ones(d)
Name = ['mzNormal','nzNormal','unNormal','mixNormal']

def b_diff(x):
    a = np.exp(x)
    return(a/(1+a))

def b_diff2(x):
    a = np.exp(x)
    b = (1 + a)**2
    return(a/b)

@nb.jit()
def GenerateData(s):
    if s=='mzNormal':
        X = np.random.multivariate_normal(mu_mz, Sigma, N)
    elif s=='nzNormal':
        X = np.random.multivariate_normal(mu_nz, Sigma, N)
    elif s=='unNormal':
        X = np.dot(np.random.multivariate_normal(mu_mz, Sigma, N), U1)
    elif s=='mixNormal':
        idx = np.random.binomial(1, p = 0.5)
        if idx == 0:
            X = np.random.multivariate_normal(mu_nz, Sigma, N)
        else:
            X = np.random.multivariate_normal(-mu_nz, Sigma, N)
    else:
        print('not founded')
    pr = (1 - 1/(1+np.exp(X.dot(beta0))))
    Y = np.zeros(N)
    for i in range(N):
        Y[i] = np.random.binomial(1, p = pr[i])
    data = np.insert(X, 0, values = Y, axis = 1)
    return(data)


def getGradient(X, y, beta, pi, r, method):
    if method == 'weighted':
        a = (b_diff(X.dot(beta)) - y)/(N * pi)
        psi = np.dot(X.T, a)/r
        return(psi)
    elif method == 'unweighted':
        a = (b_diff(X.dot(beta)) - y)
        psi = np.dot(X.T, a)/r
        return(psi)
    else:
        print('Method not found')
        

def getNewtonMatrix(X, y, beta, pi, r, method):
    if method == 'weighted':
        B = b_diff2(np.dot(X, beta))/(N * pi)
        Phi = (1/r) * (np.dot(X.T*B, X))
        return(Phi)
    elif method == 'unweighted':
        B = b_diff2(np.dot(X, beta))
        Gamma = (1/r) * (np.dot(X.T*B, X))
        return(Gamma)
    else:
        print('Method not found')
        
#@nb.jit()
def NewtonMethod(X, y, beta_int, pi, r, method, max_iter_count=100):
    if method == 'weighted':
        beta = beta_int
        for loop in range(max_iter_count):
            H = getNewtonMatrix(X, y, beta, pi, r, 'weighted')
            g = getGradient(X, y, beta, pi, r, 'weighted')
            try:
                beta_new = beta - solve(H, g, assume_a = 'sym')
            except ValueError:
                pass
            else:
                tlr = np.linalg.norm(beta_new - beta)
                beta = beta_new
                if tlr < 1e-8:
                    return(beta)
    elif method == 'unweighted':
        beta = beta_int
        for loop in range(max_iter_count):
            H = getNewtonMatrix(X, y, beta, pi, r, 'unweighted')
            g = getGradient(X, y, beta, pi, r, 'unweighted')
            try:
                beta_new = beta - solve(H, g, assume_a = 'sym')
            except ValueError:
                pass
            else:
                tlr = np.linalg.norm(beta_new - beta)
                beta = beta_new
                if tlr < 1e-8:
                    return(beta)
    else:
        print('Method not found')

def Psi_w(beta):
    XStar = SubSample[:,1:d+1]
    yStar = SubSample[:,0]
    a = (b_diff(XStar.dot(beta)) - yStar)/(N * PiSample)
    psi = XStar.T.dot(a)
    return(psi/r)
def Psi_uw(beta):
    XStar = SubSample[:,1:d+1]
    yStar = SubSample[:,0]
    a = b_diff(XStar.dot(beta)) - yStar
    psi = XStar.T.dot(a)
    return(psi/r)


@nb.jit()
def PilotEstimate(data,s):
    if s=='full':
        pilot = data
        X_pilot = pilot[:,1:d+1]
        y_pilot = pilot[:,0]
        
        reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
        reg.fit(pilot[:,1:d+1], pilot[:,0])
        beta1 = reg.coef_[0]
        
        B = b_diff2(np.dot(X_pilot, beta1))                           #A-optimality needed
        #Phi1 = (1/N) * (np.dot(X_pilot.T*B, X_pilot))                #A-optimality needed
        #Lambda1 = np.dot(np.dot(X_pilot.T, B), X_pilot)              #estimated MSE
        #InvPhi1 = solve(Phi1, np.eye(d), assume_a = 'sym')           #A-optimality needed
        #A = X_pilot.dot(InvPhi1)                                     #A-optimality needed
        A = X_pilot                                                   #L-optimiality needed
        NormsA = np.linalg.norm(A, ord = 2, axis = 1)
        b = np.sqrt(B)
        pi1 = b * NormsA
        psi_n = np.sum(pi1)
        pi = pi1/(psi_n)
        #Pilot = [beta1, pi, r0*Phi1, Lambda1]
        #Pilot = [beta1, pi, N*Phi1]
        Pilot = [beta1, pi]
        return(Pilot)
    elif s=='simple':
        r0=500
        index = np.random.choice(np.array(range(0,N)), r0, replace = False)
        pilot = data[index,]
        X_pilot = pilot[:,1:d+1]
        y_pilot = pilot[:,0]
        
        reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
        reg.fit(pilot[:,1:d+1], pilot[:,0])
        beta1 = reg.coef_[0]
        #beta1 = NewtonMethod(X_pilot, y_pilot, 0.5*np.zeros(d), (1/r0)*np.ones(r0), r0, 'unweighted')
        
        B = b_diff2(np.dot(X_pilot, beta1))                         #A-optimality needed
        Phi1 = (1/r0) * (np.dot(X_pilot.T*B, X_pilot))              #A-optimality needed
        #Lambda1 = np.dot(np.dot(X_pilot.T, B), X_pilot)             #estimated MSE needed
        X = data[:,1:d+1]
        #InvPhi1 = solve(Phi1, np.eye(d), assume_a = 'sym')          #A-optimality needed
        #A = X.dot(InvPhi1)                                          #A-optimality needed
        A = X                                                        #L-optimality needed
        NormsA = np.linalg.norm(A, ord = 2, axis = 1)
        b = np.sqrt(b_diff2(np.dot(X, beta1)))
        pi1 = b * NormsA
        psi_n = np.sum(pi1)
        pi = pi1/(psi_n)
        #Pilot = [beta1, pi, r0*Phi1, Lambda1]
        #Pilot = [beta1, pi, r0*Phi1]
        Pilot = [beta1, pi]
        return(Pilot)
    elif s=='case':
        y = data[:,0]
        n1 = np.sum(y)
        n0 = N-n1
        pi0 = (1/(2*n0))*np.ones(N)
        pi0[y==1] = 1/(2*n1)
        #print(n0)
        
        
        r0 = 500
        #index = np.random.choice(np.array(range(0,N)), r0, replace = False)
        index = np.random.choice(np.array(range(0,N)), size = r0, p = pi0, replace = True)
        pilot = data[index,]
        #pilot = data
        X_pilot = pilot[:,1:d+1]
        y_pilot = pilot[:,0]
        pi0Sample = pi0[index]
    
        reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
        reg.fit(pilot[:,1:d+1], pilot[:,0])
        beta1 = reg.coef_[0]
        #reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
        #reg.fit(data[:,1:d+1], data[:,0])
        #beta1 = reg.coef_[0]
        v = np.zeros(d)
        v[0] = 1
        beta1 = beta1 + np.log(n1/n0)*v
        
        #B = b_diff2(np.dot(X_pilot, beta1))/pi0Sample                #A-optimality needed
        #Phi1 = (1/(r0*N)) * (np.dot(X_pilot.T*B, X_pilot))           #A-optimality needed
        #Phi1 = (1/N) * (np.dot(np.dot(X_pilot.T, B), X_pilot))
        #Phi1 = (1/N) * (np.dot(np.dot(data[:,1:1+d].T, B), data[:,1:1+d]))
        #A1 = (y_pilot - b_diff(X_pilot.dot(beta1)))**2
        #Lambda1 = np.dot(np.dot(X_pilot.T, B), X_pilot)              #estimated MSE needed
        #Lambda1 = np.dot(np.dot(data[:,1:1+d].T, B), data[:,1:1+d])
    
        X = data[:,1:d+1]
        #InvPhi1 = solve(Phi1, np.eye(d), assume_a = 'sym')           #A-optimality needed
        #A = X.dot(InvPhi1)                                           #A-optimality needed
        #A = data[:,1:1+d].dot(InvPhi1)
        A = X                                                         #L-optimality needed
        NormsA = np.linalg.norm(A, ord = 2, axis = 1)
        b = np.sqrt(b_diff2(np.dot(X, beta1)))
        pi1 = b * NormsA
        psi_n = np.sum(pi1)
        pi = pi1/(psi_n)
        #Pilot = [beta1, pi, r0*Phi1, Lambda1]
        #Pilot = [beta1, pi, r0*Phi1]
        Pilot = [beta1, pi]
        return(Pilot)
    else:
        print('not found')
        
        
        
    

@nb.jit()
def MSEEstimate(X, y, beta, pi, r, method):
    B = b_diff2(np.dot(X, beta))
    #A = (y - b_diff(X.dot(beta)))**2
    if method == 'weighted':
        Phi = np.dot(np.dot(X.T, np.diag(B/(pi))), X)
        Lambda = np.dot(np.dot(X.T, np.diag(B/(pi * pi))), X) + r * np.dot(np.dot(X.T, np.diag(B/(pi))), X)
        return([Phi, Lambda])
    elif method == 'unweighted':
        Gamma = np.dot(np.dot(X.T, np.diag(B)), X)
        Gamma2 = np.dot(np.dot(X.T, np.diag(B)), X) + r * np.dot(np.dot(X.T, np.diag(B * pi)), X)
        return([Gamma, Gamma2])
    else:
        print('Method not found')
        

S = 500
R = np.ones(6)
MSE_w = np.ones(6)
MSE_uw = np.ones(6)
MSE_uni = np.ones(6)
#trV_w = np.ones(6)
#trV_uw = np.ones(6)
#RE_MSE = np.ones(6)
time_w = np.ones(6)
time_uw = np.ones(6)
time_uni = np.ones(6)
#pilot_methods=['full','simple','case']
pilot_methods=['simple']
#p_method=['A-opt','L-opt']

for s in pilot_methods:
    print(s)
    for k in range(4):
        #k = 2
        print(Name[k])
        for j in range(6):
            r = 1000 + 200 * j
            R[j] = r
            
            MSE_uni_one_r = np.zeros(S)
            time_one_uni = np.zeros(S)
            for i in range(S):
                data = GenerateData(Name[k])
                start_uni = time.time()
                index2 = np.random.choice(np.array(range(0,N)), size = r, replace = True)
                SubSample = data[index2,]
                reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
                reg.fit(SubSample[:,1:d+1], SubSample[:,0])
                beta_uni = reg.coef_[0]
                #beta_uw = NewtonMethod(SubSample[:,1:d+1], SubSample[:,0], beta1, PiSample, r, 'unweighted')
                #beta_uni = NewtonMethod(SubSample[:,1:d+1], SubSample[:,0], beta1, (1/N)*np.ones(N), r, 'unweighted')
                end_uni = time.time()
                MSE_uni_one_r[i] = np.linalg.norm(beta_uni - beta0)**2
                time_one_uni[i] = end_uni - start_uni
            
            MSE_uni[j] = np.mean(MSE_uni_one_r)
            time_uni[j] = np.sum(time_one_uni)
                
                
            MSE_w_one_r = np.zeros(S)
            time_one_w = np.zeros(S)
            for i in range(S):
                data = GenerateData(Name[k])
                start_w = time.time()
                #beta1, pi, Phi1, Lambda1= PilotEstimate(data)
                beta1, pi = PilotEstimate(data, s)
                index2 = np.random.choice(np.array(range(0,N)), size = r, p = pi, replace = True)
                SubSample = data[index2,]
                PiSample = pi[index2]
                #beta_w = fsolve(Psi_w, beta1)
                #beta_w = NewtonMethod(SubSample[:,1:d+1], SubSample[:,0], beta1, PiSample, r, 'weighted')
                reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
                reg.fit(SubSample[:,1:d+1], SubSample[:,0], sample_weight=1/PiSample)
                beta_w = reg.coef_[0]
                end_w = time.time()
                #Phi, Lambda = MSEEstimate(SubSample[:,1:d+1], SubSample[:,0], beta_w, PiSample, r, 'weighted')
                #beta_w_combine = solve(Phi1 + Phi, np.dot(Phi1, beta1) + np.dot(Phi, beta_w), assume_a = 'sym')
                #MSE_w_one_r[i] = np.linalg.norm(beta_w_combine - beta0)**2
                #trV_w_one_r[i] = np.dot(np.dot(np.linalg.inv(Phi1 + Phi), Lambda1 + Lambda),np.linalg.inv(Phi1 + Phi)).trace()
                MSE_w_one_r[i] = np.linalg.norm(beta_w - beta0)**2
                #trV_w_one_r[i] = np.dot(np.dot(np.linalg.inv(Phi), Lambda),np.linalg.inv(Phi)).trace()
                time_one_w[i] = end_w - start_w
                
            MSE_w[j] = np.mean(MSE_w_one_r)
            time_w[j] = np.sum(time_one_w)
                
                
            MSE_uw_one_r = np.zeros(S)
            time_one_uw = np.zeros(S)
            for i in range(S):
                data = GenerateData(Name[k])
                start_uw = time.time()
                #beta1, pi, Phi1, Lambda1 = PilotEstimate(data,s)
                #beta1, pi, Phi1 = PilotEstimate(data,s)
                beta1, pi = PilotEstimate(data,s)
                index2 = np.random.choice(np.array(range(0,N)), size = r, p = pi, replace = True)
                SubSample = data[index2,]
                PiSample = pi[index2]
                #beta_uw = fsolve(Psi_uw, beta1)
                reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
                reg.fit(SubSample[:,1:d+1], SubSample[:,0])
                beta_uw = reg.coef_[0]
                #beta_uw = NewtonMethod(SubSample[:,1:d+1], SubSample[:,0], beta1, PiSample, r, 'unweighted')
                end_uw = time.time()
                #Gamma, Gamma2 = MSEEstimate(SubSample[:,1:d+1], SubSample[:,0], beta_uw, PiSample, r, 'unweighted')
                #beta_uw_combine = solve(Phi1 + Gamma, np.dot(Phi1, beta1) + np.dot(Gamma, beta_uw), assume_a = 'sym')
                #MSE_uw_one_r[i] = np.linalg.norm(beta_uw_combine - beta0)**2
                #trV_uw_one_r[i] = np.dot(np.dot(np.linalg.inv(Phi1 + Gamma), Lambda1 + Gamma2),np.linalg.inv(Phi1 + Gamma)).trace()
                MSE_uw_one_r[i] = np.linalg.norm(beta_uw - beta0)**2
                #trV_uw_one_r[i] = np.dot(np.dot(np.linalg.inv(Gamma), Gamma2),np.linalg.inv(Gamma)).trace()
                time_one_uw[i] = end_uw - start_uw
        
            MSE_uw[j] = np.mean(MSE_uw_one_r)
            time_uw[j] = np.sum(time_one_uw)
            
            print('r = ',R[j], ' MSE_uni = ', MSE_uni[j], 'time_uni=', time_uni[j], ' MSE_w = ', MSE_w[j], 'time_w=', time_w[j], 
                  ' MSE_uw = ', MSE_uw[j], 'time_uw=', time_uw[j], '\n')
            
            
        MSE_full_one_r = np.zeros(S)
        time_one_full = np.zeros(S)
        for i in range(S):
            data = GenerateData(Name[k])
            start_full = time.time()
            reg = LogisticRegression(fit_intercept = False, solver = 'newton-cg')
            reg.fit(data[:,1:d+1], data[:,0])
            beta_full = reg.coef_[0]
            #beta_uw = NewtonMethod(SubSample[:,1:d+1], SubSample[:,0], beta1, PiSample, r, 'unweighted')
            end_full = time.time()
            #Gamma, Gamma2 = MSEEstimate(SubSample[:,1:d+1], SubSample[:,0], beta_uw, PiSample, r, 'unweighted')
            #beta_uw_combine = solve(Phi1 + Gamma, np.dot(Phi1, beta1) + np.dot(Gamma, beta_uw), assume_a = 'sym')
            #MSE_uw_one_r[i] = np.linalg.norm(beta_uw_combine - beta0)**2
            #trV_uw_one_r[i] = np.dot(np.dot(np.linalg.inv(Phi1 + Gamma), Lambda1 + Gamma2),np.linalg.inv(Phi1 + Gamma)).trace()
            MSE_full_one_r[i] = np.linalg.norm(beta_full - beta0)**2
            #trV_uw_one_r[i] = np.dot(np.dot(np.linalg.inv(Gamma), Gamma2),np.linalg.inv(Gamma)).trace()
            time_one_full[i] = end_full - start_full
        
        MSE_full = np.mean(MSE_full_one_r)
        time_full = np.sum(time_one_full)
    
        DF = pd.DataFrame({'r':R, 'time_uni':time_uni, 'time_w':time_w, 'time_uw':time_uw, 'time_full':time_full })
        DF.to_excel('../Results/'+Name[k]+'_uncondition_ver3.20_time_Lop'+'.xls')
        print(DF)



        plt.figure(1)
        #plt.plot(R, RE_MSE, '-o', label='RelativeMSE')
        #plt.plot(R, RE_trV, '--^', label='trV')
        #plt.plot(R, np.ones(6), '--')
        #plt.xlabel('r')
        #plt.ylabel('Relative Efficiency')
        #plt.title(Name[k])
        #plt.legend()
        #plt.figure(2)
        plt.plot(R, MSE_uni, '-s', label='MSE-uniform')
        plt.plot(R, MSE_w, '-o', label='MSE-weighted')
        plt.plot(R, MSE_uw, '--^', label = 'MSE-unweighted')
        #plt.plot(R, trV_w, '-s', label='trV-weighted')
        #plt.plot(R, trV_uw, '--D', label='trV-unweighted')
        plt.xlabel('r')
        plt.ylabel('MSE')
        plt.title(Name[k])
        plt.legend()
        plt.show()



