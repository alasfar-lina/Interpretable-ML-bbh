import numpy as np
import yaml
filename ="./data.yaml"
stream = open(filename, 'r')
data = yaml.safe_load(stream)

def XS(kq,f,W):
    R=data[f][W]['XSA0']+data[f][W]['XSA1']*kq+data[f][W]['XSA2']*kq**2
    return R*data['kl'][W]['XSSM']

def BR(kq,f):
    Rga= data[f]['Width']['gagaA0']+ data[f]['Width']['gagaA1']*kq+data[f]['Width']['gagaA2']*kq**2
    Rtot= data[f]['Width']['totA0']+ data[f]['Width']['totA1']*kq+data[f]['Width']['totA2']*kq**2
    HiggsFullwidth = data['Higgs']['width']*Rtot
    return 2*(data['Higgs']['BRbbSM']*Rga*data['Higgs']['BRgagaSM'])/Rtot


epsggF_1btag = 0.2395
epsggF_2btag = 0.0815

epsqq_1btg=0.2373
epsqq_2btg=0.0761


def weight(kq,f,W,btag):
    eps=0
    if btag==1:
        if kq==1:
            eps=epsggF_1btag
        else:
             eps= epsqq_1btg
    if btag==2:
        if kq==1:
            eps=epsggF_2btag
        else:
             eps= epsqq_2btg
    s =  XS(kq,f,W) *BR(kq,f)*eps
    return s

print( weight(1,'kl','14TeV',1))
