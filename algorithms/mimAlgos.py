
import numpy as np
import math 
from decimal import Decimal, Context
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
import os
import query
import datetime
import time
import random
from collections import OrderedDict
import entropy


###helper functions

def noise_down(lap_noise, eps_old, eps_new):
    #assert eps_new > eps_old

    pdf = [eps_old / eps_new * np.exp((eps_old - eps_new) * abs(lap_noise)),
           (eps_new - eps_old) / (2.0 * eps_new),
           (eps_old + eps_new) / (2.0 * eps_new) * (1.0 - np.exp((eps_old - eps_new) * abs(lap_noise)))]

    p = np.random.random_sample()
    if p <= pdf[0]:
        z = lap_noise

    elif p <= pdf[0] + pdf[1]:
        z = np.log(p) / (eps_old + eps_new)

    elif p <= pdf[0] + pdf[1] + pdf[2]:
        z = np.log(p * (np.exp(abs(lap_noise) * (eps_old - eps_new)) - 1.0) + 1.0) / (eps_old - eps_new)

    else:
        z = abs(lap_noise) - np.log(1.0 - p) / (eps_new + eps_old)

    return z

def integrandF1(o,ep_1,ep_2):
    F1 = lambda z:(ep_1*ep_2/4.0) * (math.exp(ep_2*(o-z))/(ep_1+ep_2) +
                                     (math.exp(ep_2*(o-z))-math.exp(ep_1*(o-z)))/(ep_1-ep_2) +
                                     math.exp(ep_1*(o-z))/(ep_1+ep_2))
        
    #F1 = lambda z:(ep_1*ep_2/4.0) * (math.exp(ep_2*(o-z))) * ((math.exp(ep_2)+math.exp(ep_1))/(ep_2+ep_1)  +
     #                                                         (math.exp(ep_2)-math.exp(ep_1))/(ep_2-ep_1)) 
    return F1
def integrandF2(o,ep_1,ep_2):
    #F2 = lambda z:(ep_1*ep_2/4.0) * (math.exp(ep_2*(z-o))) * ((math.exp(ep_2)+math.exp(ep_1))/(ep_2+ep_1)  +
     #                                                         (math.exp(ep_2)-math.exp(ep_1))/(ep_2-ep_1))
    
    F2 = lambda z:(ep_1*ep_2/4.0) * (math.exp(ep_2*(z-o))/(ep_1+ep_2) +
                                     (math.exp(ep_2*(z-o))-math.exp(ep_1*(z-o)))/(ep_1-ep_2)*1.0 +
                                     math.exp(ep_1*(z-o))/(ep_1+ep_2))
    return F2
def probx_gt_t(o,t,ep_1,ep_2):   # o is observered noisy count with ep_1, calculates prob that observed noisy count with ep_2 will be > t
    F1=integrandF1(o,ep_1,ep_2)
    F2=integrandF2(o,ep_1,ep_2)
    if(t>=o):
        ans,er=integrate.quad(F1, t, np.inf)
    else:
        ans1,er1=integrate.quad(F2, t,o)
        ans2,er2=integrate.quad(F1, o,np.inf)
        #print(ans1,er1,ans2,er2)
        ans=ans1+ans2
    
    return ans


###Threshold Shift Algorithm
def TSLM(countsByPartition, thresholds, uncertainRegion, failingRate, epsilonMax):
    ep = np.log( 0.5 / failingRate) /uncertainRegion*1.0
    
    selected = []
    epList = []
    uncertainRegionList = []
    if(ep<epsilonMax):
        lap_noises = np.random.laplace(0, 1.0/ep, len(countsByPartition))
        for i in range(len(countsByPartition)):
            epList.append(ep)
            uncertainRegionList.append(uncertainRegion)
            if(countsByPartition[i]+lap_noises[i]>thresholds[i]-uncertainRegion):
                selected.append(i)   
    return selected,epList


###Progressive Predicate-wise Laplace Mechanism
def PPWLM(maxSteps, countsByPartition, thresholds, uncertainRegion, failingRate, epsilonMax,epStart):
    ep_est = np.log( 0.5 / (failingRate/maxSteps)) / uncertainRegion*1.0
    selected = []
    epList = np.zeros(len(countsByPartition))  #eps used by each count
    positives=[]
    diff=ep_est-epStart
    if(ep_est<epsilonMax):
        ep=epStart
        base=pow(1.0*ep_est/epStart,(1.0/(maxSteps-1)))
        for j in range(1,maxSteps):
            if(len(selected)!=0):
                prev_ep=ep
                ep=epStart*pow(base,j)
                lap_noises=[noise_down(lap_noise, prev_ep, ep) for lap_noise in lap_noises]
                shiftedUncertainRegion = np.log( (1/(2.0*failingRate/maxSteps))) / ep*1.0
                for i in range(len(countsByPartition)):
                    if(i in selected):
                        if(countsByPartition[i]+lap_noises[i]>thresholds[i]+shiftedUncertainRegion):
                            epList[i]=ep
                            positives.append(i)
                            selected.remove(i)
                        elif(countsByPartition[i]+lap_noises[i]>thresholds[i]-shiftedUncertainRegion):
                            epList[i]=ep 
                        else:
                            selected.remove(i)
                            epList[i]=ep 
    return positives+selected, epList



def chooseEpsfailingRateEstimatedEntropy(j,ep_start,maxSteps,maxFineSteps,failingRate,remainingfailingRate,selected,epList,ep_prev,ep_est,lap_noises,countsByPartition,thresholds):
    beta_counter=1
    failingRate_j=failingRate/(maxSteps*1.0)
    remainingfailingRate=round(remainingfailingRate,6)
    base=pow(1.0*ep_est/ep_start,(1.0/(maxSteps-1)))
    ep_j=ep_start*pow(base,j)
    ep_jminus1=ep_prev
    ep=ep_prev
    opt_j=j
    opt_ent=0
    ep_opt=ep_est
    optfailingRate = remainingfailingRate       
    while(ep<=ep_est and remainingfailingRate>0 and remainingfailingRate>((failingRate/maxSteps))): 
        for fs in range(maxFineSteps):
            ep_range=ep_j-ep_jminus1
            if(ep_range<=0):
                break
            ep=ep+1.0*ep_range/maxFineSteps
            exp_pos=0
            for i in range(len(lap_noises)):
                if(i in selected):
                    shiftedUncertainRegion = np.log( 0.5 / failingRate_j) / ep*1.0
                    prob_pos_eliminated=probx_gt_t(countsByPartition[i]+lap_noises[i],thresholds[i]+shiftedUncertainRegion,ep_prev,ep)
                    prob_pos_selected=probx_gt_t(countsByPartition[i]+lap_noises[i],thresholds[i]-shiftedUncertainRegion,ep_prev,ep)
                    exp_pos+=(prob_pos_selected-prob_pos_eliminated)
            cost=list(epList)
            cost=sorted(cost)
            for i in selected:
                cost.remove(ep_prev)
                cost.append(ep)
            if math.isnan(exp_pos):
                exp_pos=0
            for i in range(0,int(round(exp_pos))):
                cost.remove(ep)
                cost.append(ep_est)
            found=dict_ent.get(str(cost))
            if   found==-1:  
                new_ent=entropy.computeMinEntropy(cost)[0]
                dict_ent.put(str(cost),new_ent)
                #print(new_ent)
            else:
                matches[0]+=1
                new_ent=found
            
            if(opt_ent<new_ent and int(round(exp_pos))!=len(selected)):
                opt_ent=new_ent
                ep_opt=ep
                optfailingRate=1.0*(failingRate/maxSteps)*beta_counter
                opt_j=j
        beta_counter+=1
        j+=1
        ep_jminus1=ep_j
        ep_j=ep_start*pow(base,j)   
    return ep_opt, round(optfailingRate,6),opt_j

###Data Dependent Predicate-wise Laplace Mechanism
def DDPWLM(maxSteps, maxFineSteps,ep_start, countsByPartition, thresholds, uncertainRegion, failingRate, epsilonMax): 
    ep_est = np.log( 0.5 / (failingRate/maxSteps)) / uncertainRegion*1.0
    remainingfailingRate=failingRate
    selected=[]
    positives=[]
    epList = np.zeros(len(countsByPartition))  #eps used by each count
    if(ep_est<epsilonMax):
        ep = ep_start
        lap_noises = np.random.laplace(0, 1.0/ep, len(countsByPartition))
        for i in range(len(countsByPartition)):
            shiftedUncertainRegion = np.log( 0.5*maxSteps / failingRate) / ep*1.0
            if(countsByPartition[i]+lap_noises[i]>thresholds[i]+shiftedUncertainRegion):
                positives.append(i)
            elif(countsByPartition[i]+lap_noises[i]>thresholds[i]-shiftedUncertainRegion):
                selected.append(i)
            epList[i] = ep  
        j=1
        remainingfailingRate-=failingRate/maxSteps*1.0
        while (ep<ep_est and round(remainingfailingRate,6)>0 and len(selected)!=0):
            j=j+1
            prev_ep=ep
            ep,failingRate_j,j_traversed=(j,ep_start,maxSteps,maxFineSteps,failingRate,remainingfailingRate,selected,epList,ep_prev,ep_est,lap_noises,countsByPartition,thresholds)
            lap_noises=[noise_down(lap_noise, prev_ep, ep) for lap_noise in lap_noises]
            for i in range(len(countsByPartition)):
                shiftedUncertainRegion = np.log( 0.5 / failingRate_j) / ep*1.0
                if(i in selected):
                    if( countsByPartition[i]+lap_noises[i]>thresholds[i]+shiftedUncertainRegion):
                        epList[i]=ep 
                        positives.append(i)
                        selected.remove(i)
                    elif( countsByPartition[i]+lap_noises[i]>thresholds[i]-shiftedUncertainRegion):
                        epList[i]=ep 
                    else:
                        selected.remove(i)
                        epList[i]=ep
            remainingfailingRate-= failingRate_j
            remainingfailingRate = round(remainingfailingRate,6)
    return positives+selected, epList 

