import collections
import operator
import math
from bisect import bisect_left
import sys
import os
import glob
import commonfn

param_set={
        "mzML_files",
        "min_highest_I",
        "length_of_ion_chromatogram",
        }
param_dict=commonfn.read_param(param_set)

mzML_files=glob.glob(param_dict["mzML_files"])


min_feat_height=float(param_dict["min_highest_I"])
min_auc=min_feat_height*2 #float(param_dict["min_auc"])
peak_w=[int(x) for x in param_dict["length_of_ion_chromatogram"].split()]

Point=collections.namedtuple('Point',('rt mz I'))
Coef=collections.namedtuple('Coef',('rt sc coef'))
Peak=collections.namedtuple('Peak',('mz rt sc coef auc'))

wave_scales=[x/2 for x in list(range(peak_w[1]+2,peak_w[0]-3,-2))]
wave_scales=[1.5]+list(range(2,25))+list(range(25,99,2))
ws0=bisect_left(wave_scales,peak_w[0]/2)
ws1=bisect_left(wave_scales,peak_w[1]/2,lo=ws0)
wave_scales=wave_scales[:ws1+1+1][::-1]
wave_sqrt=[math.sqrt(x) for x in wave_scales]

def cwt(mzML_file):
    basename0=os.path.basename(mzML_file)

    filei='eic_'+basename0+'.txt'

    def get_EICs():
        eic_dat=open(filei)
        for line in eic_dat:
            lsp=line.rstrip('\n').split('\t')
            if len(lsp)==3:
                EIC.append(Point(*[float(x) for x in lsp]))
            elif line=='-\n':
                yield (EIC,rt_all)
                EIC=[]
            elif line.startswith('scan '):
                EIC=[]
                rt_all=sorted(float(x) for x in next(eic_dat).rstrip().split())

    def findridge(eic):
        EIC,rt_all=eic
        eic_dict={pt.rt:(pt.mz,pt.I) for pt in EIC}
        eic_rt=sorted(x for x,y in eic_dict.items() if y[1]>min_feat_height)#reduce exe time

        coefs = [[0]*len(eic_rt) for i in wave_scales]
        for xx,wave_scale in enumerate(wave_scales):
            for yy,wave_location in enumerate(eic_rt):
                bd=max(wave_scale,9)
                pos0=bisect_left(rt_all,wave_location-bd)
                pos1=bisect_left(rt_all,wave_location+bd,lo=pos0)
                rt_=rt_all[pos0:pos1]
                int_I=[0]*len(rt_)
                for i,rt0 in enumerate(rt_):
                    if rt0 in eic_dict:
                        tsig2=((rt0 - wave_location)/wave_scale)**2
                        int_I[i]=eic_dict[rt0][1]*math.exp(-tsig2/2)*(1.-tsig2)
                coefs[xx][yy]=sum((I0+I1)*(rt1-rt0) for rt0,rt1,I0,I1 in zip(rt_,rt_[1:],int_I,int_I[1:]))/2/wave_sqrt[xx]


        rtdict={rt_i:n for n,rt_i in enumerate(eic_rt)}
        max_map = [[0]*len(eic_rt) for i in wave_scales]
        for xx,wave_scale in enumerate(wave_scales):
            max_coefs=[]
            coef_xx=[(eic_rt[n],x) for n,x in enumerate(coefs[xx]) if x>0]
            while any(y>0 for _,y in coef_xx):
                max_coef=max(coef_xx,key=operator.itemgetter(1))
                max_i=coef_xx.index(max_coef)
                if (max_i!=0 and max_i!=len(coef_xx)-1) and (coef_xx[max_i-1][1]==0 or coef_xx[max_i+1][1]==0):
                    coef_xx[max_i]=(coef_xx[max_i][0],0)
                    continue
                max_coefs.append(max_coef)
                for i in range(len(coef_xx)):
                    if abs(max_coef[0]-coef_xx[i][0])<max(wave_scale*1.5,8):
                        coef_xx[i]=(coef_xx[i][0],0)
            for rt,coef in max_coefs:
                max_map[xx][rtdict[rt]]=coef


        ridgelines=[]
        for xx,scale_coef in enumerate(max_map):
            for yy,coef in enumerate(scale_coef):
                added=False
                if coef:
                    for rl in ridgelines:
                        if ((abs(bisect_left(rt_all,rl[-1].rt)-bisect_left(rt_all,eic_rt[yy]))<2 or abs(rl[-1].rt-eic_rt[yy])<1) and \
                                rl[-1].sc==wave_scales[xx-1]):
                            rl.append(Coef(eic_rt[yy],wave_scales[xx],coef))
                            added=True
                            break
                    if not added:
                        ridgelines.append([Coef(eic_rt[yy],wave_scales[xx],coef)])

        ridgelines=[x for x in ridgelines if len(x)>8]




        for rd in ridgelines[:]:
            peak_loc=max(rd,key=operator.attrgetter('coef'))
            if not peak_w[0]<=peak_loc.sc*2<=peak_w[1]:
                ridgelines.remove(rd)
                continue

            pos0=bisect_left(rt_all,peak_loc.rt-peak_loc.sc)
            pos1=bisect_left(rt_all,peak_loc.rt+peak_loc.sc,lo=pos0)
            if pos0==0 or pos1==len(rt_all) or sum(1 for rt0 in rt_all[pos0:pos1] if rt0 in eic_dict)/(pos1-pos0)<.5:
                ridgelines.remove(rd)
                continue

            pos=bisect_left(rt_all,peak_loc.rt,lo=pos0,hi=pos1)
            apex=[eic_dict[rt0][1] for rt0 in rt_all[pos-2:pos+3] if rt0 in eic_dict]
            left0=[eic_dict[rt0][1] for rt0 in rt_all[pos0-2:pos] if rt0 in eic_dict]
            right0=[eic_dict[rt0][1] for rt0 in rt_all[pos+1:pos1+2] if rt0 in eic_dict]
            if len(apex)>1 and left0 and right0:
                apex_h=sorted(apex)[-2]/2
                if apex_h<min(left0) or apex_h<min(right0):
                    ridgelines.remove(rd)
                    continue
            else:
                ridgelines.remove(rd)
                continue

        peaks=[]

        for rd in ridgelines:
            peak_loc=max(rd,key=operator.attrgetter('coef'))
            pos0=bisect_left(rt_all,peak_loc.rt-1.5*peak_loc.sc)
            pos1=bisect_left(rt_all,peak_loc.rt+1.5*peak_loc.sc,lo=pos0)
            auc=sum((eic_dict.get(rt0,(0,0))[1]+eic_dict.get(rt1,(0,0))[1])*(rt1-rt0) for rt0,rt1 in zip(rt_all[pos0:],rt_all[pos0+1:pos1]))/2
            if min_auc<auc:# and peak_loc.sc<=len(rd):#ridge length
                rt_sub=[rt for rt in rt_all[pos0:pos1] if rt in eic_dict]
                peakmz=sum(eic_dict[rt][0]*eic_dict[rt][1] for rt in rt_sub)/sum(eic_dict[rt][1] for rt in rt_sub)
                peaks.append(Peak(peakmz,peak_loc.rt,peak_loc.sc,peak_loc.coef,auc))
        return peaks


    with open('ms1feature_'+basename0+'.txt','w') as writepeak:

        peak_list=[]
        for peaks in map(findridge, get_EICs()):
            peak_list.extend(peaks)

        peak_list.sort()
        for peak in peak_list[:]:
            if peak in peak_list:
                peak_mz=sorted((x for x in peak_list[bisect_left(peak_list,(peak.mz,)):bisect_left(peak_list,(peak.mz+.01,))] if abs(x.rt-peak.rt)<x.sc+peak.sc),key=operator.attrgetter('coef'),reverse=True)
                for peak0 in peak_mz[1:]:
                    peak_list.remove(peak0)


        writepeak.write('MS1\n')
        for peak in peak_list:
            writepeak.write('\t'.join(str(x) for x in peak)+'\n')

