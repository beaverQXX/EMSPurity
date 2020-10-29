import numpy as np 
import pandas as pd 
import matplotlib.pylab as plt
from scipy import stats
import math
import warnings 
import random
warnings.filterwarnings("ignore") 
def getDataSet():
    with open("netural_SNP.txt") as fb:
        lines = fb.readlines()
    data = []
    for line in lines[1:]:
        line = line.strip().split('\t')
        value = float(line[8])/float(line[6])
        if(value < 0.96):
            data.append(float(line[8])/float(line[6]))
    return data
def getLogLikelihood(w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,sim_datas):
    loglikelihood = 0.0
    rz1_up = w1 * stats.norm(mu1,std1).pdf(sim_datas)
    rz2_up = w2 * stats.norm(mu2, std2).pdf(sim_datas)
    rz3_up = w3 * stats.norm(mu3, std3).pdf(sim_datas)
    rz4_up = w4 * stats.norm(mu4, std4).pdf(sim_datas)
    for i in range(len(sim_datas)):
        loglikelihood += math.log(rz1_up[i]+rz2_up[i]+rz3_up[i]+rz4_up[i]+1e-5)
    return loglikelihood
def getMiu(miu1):
    miu2 = 2*miu1
    miu3 = 0.5
    miu4 = miu3+miu1
    print("miu1",miu1)
    return miu1,miu2,miu3,miu4
def EM(w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,sim_datas):
    n = len(sim_datas)  # 样本长度
    preLogLikehood = getLogLikelihood(w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,sim_datas)
    prew1 = w1;prew2 = w2;prew3 = w3;prew4 = w4;
    premu1 = mu1;premu2 = mu2;premu3 = mu3;premu4 = mu4;
    prestd1 = std1;prestd2 = std2;prestd3 = std3;prestd4 = std4;
    flag  = True
    d = 1

    # 开始EM算法的主循环
    for t in range(100):
        # E-step
        rz1_up = w1 * stats.norm(mu1,std1).pdf(sim_datas)
        rz2_up = w2 * stats.norm(mu2, std2).pdf(sim_datas)
        rz3_up = w3 * stats.norm(mu3, std3).pdf(sim_datas)
        rz4_up = w4 * stats.norm(mu4, std4).pdf(sim_datas)
        rz_down = rz1_up+rz2_up+rz3_up+rz4_up
        if(np.isnan(rz_down[0])):
            flag = False
            return prew1,prew2,prew3,prew4,premu1,premu2,premu3,premu4,prestd1,prestd2,prestd3,prestd4,preLogLikehood
       # print(rz_down)
        rz1 = rz1_up/rz_down      # 为分布1的概率
        rz2 = rz2_up/rz_down      # 为分布2的概率
        rz3 = rz3_up/rz_down      # 为分布1的概率
        rz4 = rz4_up/rz_down
        w1 = np.sum(rz1)/n
        w2 = np.sum(rz2)/n
        w3 = np.sum(rz3)/n
        w4 = np.sum(rz4)/n
            # M-step
        mu1 = np.sum(rz1*sim_datas)/np.sum(rz1)
        mu2 = np.sum(rz2*sim_datas)/np.sum(rz2)
        mu3 = np.sum(rz3*sim_datas)/np.sum(rz3)
        mu4 = np.sum(rz4*sim_datas)/np.sum(rz4)
        std1 = np.sqrt(np.sum(rz1*np.square(sim_datas-mu1))/(d*np.sum(rz1)))
        std2 = np.sqrt(np.sum(rz2*np.square(sim_datas-mu2))/(d*np.sum(rz2)))
        std3 = np.sqrt(np.sum(rz3*np.square(sim_datas-mu3))/(d*np.sum(rz3)))
        std4 = np.sqrt(np.sum(rz4*np.square(sim_datas-mu4))/(d*np.sum(rz4)))
        loglikelihood = getLogLikelihood(w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,sim_datas)
        print(t+1,":",loglikelihood)
        if(np.isnan(loglikelihood)):
            flag = False
            return prew1,prew2,prew3,prew4,premu1,premu2,premu3,premu4,prestd1,prestd2,prestd3,prestd4,preLogLikehood
        if(abs(preLogLikehood-loglikelihood)<1e-5):
            break
        else:
            preLogLikehood = loglikelihood
        prew1 = w1;prew2 = w2;prew3 = w3;prew4 = w4;
        premu1 = mu1;premu2 = mu2;premu3 = mu3;premu4 = mu4;
        prestd1 = std1;prestd2 = std2;prestd3 = std3;prestd4 = std4;
    return w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,preLogLikehood
def getMaxTwo(w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,logLikehood,
    sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_1_logLikehood,
    sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4,sample_2_logLikehood ):
    #3最小
    if(min(logLikehood,sample_1_logLikehood) >= sample_2_logLikehood):
        return w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4
    #2最小
    elif(min(logLikehood ,sample_2_logLikehood) >= sample_1_logLikehood):
        return w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4
    #1最小
    elif(min(sample_1_logLikehood ,logLikehood) >= logLikehood):
        return sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4

def getMaxTwo2(
    s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,s_1_logLikehood,
    s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4,s_2_logLikehood,
    sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_1_logLikehood,
    sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4,sample_2_logLikehood ):
    #12/34 12大或34
    if(s_1_logLikehood>=max(sample_1_logLikehood,sample_2_logLikehood) and s_2_logLikehood>=max(sample_1_logLikehood,sample_2_logLikehood)):
        return s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4
    elif(sample_1_logLikehood>=max(s_1_logLikehood,s_2_logLikehood) and sample_2_logLikehood>=max(s_1_logLikehood,s_2_logLikehood)):
        return sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4
    #13/24 13大或24大
    elif(s_1_logLikehood>=max(s_2_logLikehood,sample_2_logLikehood) and sample_1_logLikehood>=max(s_2_logLikehood,sample_2_logLikehood)):
        return s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4
    elif(s_2_logLikehood>=max(s_1_logLikehood,sample_1_logLikehood) and sample_2_logLikehood>=max(s_1_logLikehood,sample_1_logLikehood)):
        return s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4

    #14/23 14大或23大
    elif(s_1_logLikehood>=max(s_2_logLikehood,sample_1_logLikehood) and sample_2_logLikehood>=max(s_2_logLikehood,sample_1_logLikehood)):
        return s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4
    elif(s_2_logLikehood>=max(s_1_logLikehood,sample_2_logLikehood) and sample_1_logLikehood>=max(s_1_logLikehood,sample_2_logLikehood)):
        return s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4

sim_datas = getDataSet()#获得实验数据
miu1 = 0.2
mu1,mu2,mu3,mu4 =  getMiu(miu1)
candidate = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45]
index = np.random.randint(0,9,50)
randNum = 0
w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,logLikehood=EM(0.25,0.25,0.25,0.25,mu1,mu2,mu3,mu4,0.01,0.01,0.01,0.01,sim_datas)
#第一个随机搜索点
sample_1_miu1 = candidate[index[randNum]]
randNum += 1
sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4= getMiu(sample_1_miu1)
sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_1_logLikehood=EM(0.25,0.25,0.25,0.25,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,0.01,0.01,0.01,0.01,sim_datas)

#第二个随机搜索点
sample_2_miu1 = candidate[index[randNum]]
randNum+=1
sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4 =  getMiu(sample_2_miu1)
sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4,sample_2_logLikehood=EM(0.25,0.25,0.25,0.25,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,0.01,0.01,0.01,0.01,sim_datas)
print("loglikelihood:",sample_1_logLikehood,sample_2_logLikehood)
s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4=getMaxTwo(w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,logLikehood,sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_1_logLikehood,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4,sample_2_logLikehood)
res_w1=0;res_w2=0;res_w3=0;res_w4=0;res_miu1=0;res_miu2=0;res_miu3=0;res_miu4=0;res_std1=0;res_std2=0;res_std3=0;res_std4=0
while(randNum<50):
    #进行EM算法
    s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,s_1_logLikehood=EM(s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,sim_datas)
    #进行EM算法
    s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4,s_2_logLikehood=EM(s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4,sim_datas)

    sample_1_miu1 = candidate[index[randNum]]
    randNum += 1
    sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4= getMiu(sample_1_miu1)
    sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_1_logLikehood=EM(0.25,0.25,0.25,0.25,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,0.01,0.01,0.01,0.01,sim_datas)
    sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_1_logLikehood=EM(0.25,0.25,0.25,0.25,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,0.01,0.01,0.01,0.01,sim_datas)

    #第二个随机搜索点
    sample_2_miu1 = candidate[index[randNum]]
    randNum+=1
    sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4 =  getMiu(sample_2_miu1)
    sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4,sample_2_logLikehood=EM(0.25,0.25,0.25,0.25,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,0.01,0.01,0.01,0.01,sim_datas)
    s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4=getMaxTwo2(s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4,s_1_logLikehood,s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4,s_2_logLikehood,sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4,sample_1_logLikehood,sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4,sample_2_logLikehood )
    randNum += 1
    if(randNum == 50):
        print("randNum:",50)
        print(s_1_logLikehood,s_2_logLikehood,sample_1_logLikehood,sample_2_logLikehood)
        if(s_1_logLikehood == max(s_1_logLikehood,s_2_logLikehood,sample_1_logLikehood,sample_1_logLikehood)):
            res_w1,res_w2,res_w3,res_w4,res_miu1,res_miu2,res_miu3,res_miu4,res_std1,res_std2,res_std3,res_std4=s_1_w1,s_1_w2,s_1_w3,s_1_w4,s_1_miu1,s_1_miu2,s_1_miu3,s_1_miu4,s_1_std1,s_1_std2,s_1_std3,s_1_std4
        elif(s_2_logLikehood == max(s_1_logLikehood,s_2_logLikehood,sample_1_logLikehood,sample_2_logLikehood)):
            res_w1,res_w2,res_w3,res_w4,res_miu1,res_miu2,res_miu3,res_miu4,res_std1,res_std2,res_std3,res_std4=s_2_w1,s_2_w2,s_2_w3,s_2_w4,s_2_miu1,s_2_miu2,s_2_miu3,s_2_miu4,s_2_std1,s_2_std2,s_2_std3,s_2_std4
        elif(sample_1_logLikehood == max(s_1_logLikehood,s_2_logLikehood,sample_1_logLikehood,sample_2_logLikehood)):
            res_w1,res_w2,res_w3,res_w4,res_miu1,res_miu2,res_miu3,res_miu4,res_std1,res_std2,res_std3,res_std4=sample_1_w1,sample_1_w2,sample_1_w3,sample_1_w4,sample_1_miu1,sample_1_miu2,sample_1_miu3,sample_1_miu4,sample_1_std1,sample_1_std2,sample_1_std3,sample_1_std4
        elif(sample_2_logLikehood == max(s_1_logLikehood,s_2_logLikehood,sample_1_logLikehood,sample_2_logLikehood)):
            res_w1,res_w2,res_w3,res_w4,res_miu1,res_miu2,res_miu3,res_miu4,res_std1,res_std2,res_std3,res_std4=sample_2_w1,sample_2_w2,sample_2_w3,sample_2_w4,sample_2_miu1,sample_2_miu2,sample_2_miu3,sample_2_miu4,sample_2_std1,sample_2_std2,sample_2_std3,sample_2_std4

    #随机搜索
print(res_w1,res_w2,res_w3,res_w4)
print(res_miu1,res_miu2,res_miu3,res_miu4)
res_miu2 = 2*res_miu1
res_miu3 = 0.5
res_miu4 = res_miu3+res_miu1
w1,w2,w3,w4,mu1,mu2,mu3,mu4,std1,std2,std3,std4,loglikelihood=EM(0.25,0.25,0.25,0.25,res_miu1,res_miu2,res_miu3,res_miu4,0.01,0.01,0.01,0.01,sim_datas)
print("mu1",mu1)
print("mu2",mu2)
print("mu3",mu3)
print("mu4",mu4)
print(loglikelihood)
predict11 = stats.norm(mu1, std1).pdf(sim_datas)
predict22 = stats.norm(mu2, std2).pdf(sim_datas)
predict33 = stats.norm(mu3, std3).pdf(sim_datas)
predict44 = stats.norm(mu4, std4).pdf(sim_datas)
purity = []
for i in range(len(sim_datas)):
    # if(predict22[i]>predict11[i] and predict22[i]>predict33[i] and predict22[i]>predict44[i]):
    #     purity.append(sim_datas[i])
    # el
    if(predict11[i]>predict22[i] and predict11[i]>predict33[i] and predict11[i]>predict44[i]):
        purity.append(sim_datas[i]*2)
    elif(predict44[i]>predict11[i] and predict44[i]>predict22[i] and predict44[i]>predict33[i]):
        purity.append(sim_datas[i])
print(np.mean(purity))




