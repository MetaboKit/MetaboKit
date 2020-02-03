import collections
import operator
import math
from bisect import bisect_left
import sys
#import concurrent.futures
import os
import glob
#from statistics import median
import commonfn

param_set={
        "mzML_files",
        "min_highest_I",
        #"min_auc",
        "length_of_ion_chromatogram",
        #"num_threads",
        }
param_dict=commonfn.read_param(param_set)

mzML_files=glob.glob(param_dict["mzML_files"])

#num_threads=2#int(param_dict["num_threads"])

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
#wave_scales=wave_scales[max(ws0-1,0):ws1+2][::-1]
wave_scales=wave_scales[:ws1+1+1][::-1]
wave_sqrt=[math.sqrt(x) for x in wave_scales]
#print(min(wave_scales),max(wave_scales))
#sys.exit()
#def ricker_wavelet(x,scalParam):
#    tsig2=(x/scalParam)**2
#    #return math.exp(-tsig2/2.)/math.sqrt(scalParam)*(1.-tsig2)
#    return math.exp(-x**2 / (2.0 * scalParam**2)) * 2. / math.sqrt(3. * scalParam * math.sqrt(math.pi)) * (1. - (x/scalParam)**2)
#sc=1.1
#print(ricker_wavelet(sc*math.sqrt(3)-.0000,sc))#minima
#print(-4/math.sqrt(3*sc)/math.pow(math.pi,1/4)/math.exp(1.5))
#print(ricker_wavelet(0,sc))
#print(2/math.sqrt(3*sc)/math.pow(math.pi,1/4))
#sys.exit()
#def ricker_wave(x_scalParam):
#    tsig2=(x_scalParam)**2
#    return math.exp(-tsig2/2.)*(1.-tsig2)

#print(len(mzML_files),'mzML files')
#for mzML_file in mzML_files:
def cwt(mzML_file):
    basename0=os.path.basename(mzML_file)
    #print(basename0)

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

    #for eic in EICs[956:957]:
    def findridge(eic):
        EIC,rt_all=eic
        #eic_dict=dict() # highest intensity in this m/z range
        #for pt in EIC:
        #    if pt.rt not in eic_dict or eic_dict[pt.rt][1]<pt.I:
        #        eic_dict[pt.rt]=(pt.mz,pt.I)
        eic_dict={pt.rt:(pt.mz,pt.I) for pt in EIC}
        eic_rt=sorted(x for x,y in eic_dict.items() if y[1]>min_feat_height)#reduce exe time
        #eic_rt=sorted(eic_dict.keys())
        #rt_I=sorted(pt for pt in eic_dict.values())

        coefs = [[0]*len(eic_rt) for i in wave_scales]
        for xx,wave_scale in enumerate(wave_scales):
            for yy,wave_location in enumerate(eic_rt):
                bd=max(wave_scale,8)
                pos0=bisect_left(rt_all,wave_location-bd)
                pos1=bisect_left(rt_all,wave_location+bd,lo=pos0)
                rt_=rt_all[pos0:pos1]
                int_I=[0]*len(rt_)
                for i,rt0 in enumerate(rt_):
                    if rt0 in eic_dict:
                        tsig2=((rt0 - wave_location)/wave_scale)**2
                        int_I[i]=eic_dict[rt0][1]*math.exp(-tsig2/2)*(1.-tsig2)
                coefs[xx][yy]=sum((I0+I1)*(rt1-rt0) for rt0,rt1,I0,I1 in zip(rt_,rt_[1:],int_I,int_I[1:]))/2/wave_sqrt[xx]

        ###print coef
        #if 363.3<sum([pt.mz for pt in EIC])/len(EIC)<363.33:
        #    print(sum([pt.mz for pt in EIC])/len(EIC))
        #    for x,w in zip(coefs,wave_scales):
        #        print(w,[(round(xx,1),round(yy,1)) for xx,yy in zip(eic_rt,x) if 0<xx<2000])

        rtdict={rt_i:n for n,rt_i in enumerate(eic_rt)}
        # mark local maximum in cwt coeficient matrix
        max_map = [[0]*len(eic_rt) for i in wave_scales]
        for xx,wave_scale in enumerate(wave_scales):
            max_coefs=[]
            coef_xx=[(eic_rt[n],x) for n,x in enumerate(coefs[xx]) if x>0]
            while any(y>0 for _,y in coef_xx):
                # 2.
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
        #sys.exit()


        #construct ridgelines
        ridgelines=[]
        for xx,scale_coef in enumerate(max_map):
            for yy,coef in enumerate(scale_coef):
                added=False
                if coef:
                    for rl in ridgelines:
                        #if (abs(rl[-1].rt-eic_rt[yy])<2 and \
                        if ((abs(bisect_left(rt_all,rl[-1].rt)-bisect_left(rt_all,eic_rt[yy]))<2 or abs(rl[-1].rt-eic_rt[yy])<1) and \
                                rl[-1].sc==wave_scales[xx-1]):
                            rl.append(Coef(eic_rt[yy],wave_scales[xx],coef))
                            added=True
                            break
                    if not added:
                        ridgelines.append([Coef(eic_rt[yy],wave_scales[xx],coef)])

        #filter ridgeline for length
        ridgelines=[x for x in ridgelines if len(x)>8]
        #ridgelines=[x for x in ridgelines if x[0].sc-x[-1].sc>5/60]
        #ridgelines=[x for x in ridgelines if max(y[2] for y in x)>0]

        #if 363.3<sum([pt.mz for pt in EIC])/len(EIC)<363.33:
        #    print('len ridgelines',len(ridgelines))
        #    for rd in ridgelines: print(rd)



        ## wavelet coef filter, peak property filter
        for rd in ridgelines[:]:
            #print(len(rd))
            peak_loc=max(rd,key=operator.attrgetter('coef'))
            #if math.isclose(peak_loc.sc,wave_scales[0]) or math.isclose(peak_loc.sc,wave_scales[-1]):
            # wavelet coef filter
            if not peak_w[0]<=peak_loc.sc*2<=peak_w[1]:
                ridgelines.remove(rd)
                continue

            #scan ratio filter
            pos0=bisect_left(rt_all,peak_loc.rt-peak_loc.sc)
            pos1=bisect_left(rt_all,peak_loc.rt+peak_loc.sc,lo=pos0)
            if pos0==0 or pos1==len(rt_all) or sum(1 for rt0 in rt_all[pos0:pos1] if rt0 in eic_dict)/(pos1-pos0)<.5:
                ridgelines.remove(rd)
                continue

            # background filter
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
            ##peak property filter
            pos0=bisect_left(rt_all,peak_loc.rt-peak_loc.sc)
            pos1=bisect_left(rt_all,peak_loc.rt+peak_loc.sc,lo=pos0)
            auc=sum((eic_dict.get(rt0,(0,0))[1]+eic_dict.get(rt1,(0,0))[1])*(rt1-rt0) for rt0,rt1 in zip(rt_all[pos0:],rt_all[pos0+1:pos1]))/2
            #auc =eic_dict.get(rt_all[pos0],(0,0))[1]*(rt_all[pos0+1]-rt_all[pos0])
            #auc+=eic_dict.get(rt_all[pos1-1],(0,0))[1]*(rt_all[pos1-1]-rt_all[pos1-2])
            #for p in range(pos0+1,pos1-1):
            #    auc+=eic_dict.get(rt_all[p],(0,0))[1]*(rt_all[p+1]-rt_all[p-1])
            #auc/=2
            if min_auc<auc:# and peak_loc.sc<=len(rd):#ridge length
                rt_sub=[rt for rt in rt_all[pos0:pos1] if rt in eic_dict]
                peakmz=sum(eic_dict[rt][0]*eic_dict[rt][1] for rt in rt_sub)/sum(eic_dict[rt][1] for rt in rt_sub)
                #peaks.append(Peak(peakmz,peak_loc.rt,peak_loc.sc,peak_loc.coef/auc,auc))
                peaks.append(Peak(peakmz,peak_loc.rt,peak_loc.sc,peak_loc.coef,auc))
        return peaks


    with open('ms1feature_'+basename0+'.txt','w') as writepeak:

        peak_list=[]
        #with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        #    for peaks in executor.map(findridge, EICs):
        #        peak_list.extend(peaks)
        for peaks in map(findridge, get_EICs()):
            peak_list.extend(peaks)

        peak_list.sort()
        for peak in peak_list[:]:
            #if peak not yet removed
            #if peak_list[bisect_left(peak_list,peak)] is peak:
            if peak in peak_list:
                peak_mz=sorted((x for x in peak_list[bisect_left(peak_list,(peak.mz,)):bisect_left(peak_list,(peak.mz+.01,))] if abs(x.rt-peak.rt)<x.sc+peak.sc),key=operator.attrgetter('coef'),reverse=True)
                for peak0 in peak_mz[1:]:
                    peak_list.remove(peak0)

        #print('len peak_list',len(peak_list))
        #print(peak_list==sorted(peak_list))

        writepeak.write('MS1\n')
        for peak in peak_list:
            writepeak.write('\t'.join(str(x) for x in peak)+'\n')

