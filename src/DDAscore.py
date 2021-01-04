
import xml.etree.ElementTree as ET
import collections
import operator
import sys
import math
from bisect import bisect_left
import re
import os
import glob
import concurrent.futures

import DDAreadlib
import DDAcommonfn
from DDAcommonfn import bound_ppm

Point=collections.namedtuple('Point',('mz rt I'))
Spec=collections.namedtuple('Spec',('ms1mz rt mz I mz_all I_all'))
Peak=collections.namedtuple('Peak',('mz rt sc coef auc mmz'))
Ent=collections.namedtuple('Ent',('Mmass name mz I adduct charge rt formu'))

param_set={
        "mzML_files",
        "library",
        "adduct",
        "ms1_ppm",
        "ms2_ppm",
        "ISF_rt_diff",
        "ISF_score",
        "MS2_score",
        "min_peaks",
        "base_peak_filter",
        "RT_shift",
        }
param_dict=DDAcommonfn.read_param(param_set)

RT_shift=float(param_dict['RT_shift'])

lib_types=[x.split(' #',1)[0].strip() for x in param_dict["library"].splitlines()]

ms1ppm=float(param_dict['ms1_ppm'])/1e6
ms2ppm=float(param_dict['ms2_ppm'])/1e6

min_peaks=int(param_dict["min_peaks"])
bpfilter=param_dict.get("base_peak_filter",'0.00 20').split(' ')
bpfilter=[float(bpfilter[0]),int(bpfilter[1])]

ISF_rt_diff=float(param_dict["ISF_rt_diff"])
ISF_score=float(param_dict["ISF_score"])
MS2_score=float(param_dict["MS2_score"])

adduct_list=dict()
adduct_list['?']=(0, 1, '?')

mzML_files=sorted(glob.glob(param_dict["mzML_files"]))

for line in open(mzML_files[0]):
    if 'MS:1000130'in line:
        ispos=True
        break
    elif 'MS:1000129'in line:
        ispos=False
        break

def cos_sim(list1,list2):
    if len(list1)!=len(list2):
        print('adf')
        sys.exit()
    if sum(list1)<=0 or sum(list2)<=0: return 0
    return sum(math.sqrt(x*y) for x,y in zip(list1,list2))/math.sqrt(sum(list1)*sum(list2))
    return sum(x*y for x,y in zip(list1,list2))/math.sqrt(sum(x*x for x in list1)*sum(x*x for x in list2))


def readdatbase():
    dpath=libpaths[0]
    dlist=[]
    with open(dpath) as dbase:
        line1=dbase.readline()
        for x in (line.rstrip('\n').split('\t') for line in dbase):
            if len(lib_types[0])>9:
                dlist.append((float(x[1]),x[0]))
            else:
                dlist.append((float(x[1]),x[0]+' '+(x[2] if x[2] else x[3])))
    return sorted(dlist)

if "NoMatch" in lib_types[0]:
    dlist=readdatbase()
else:
    lib_ent=DDAreadlib.get_cpds()

def read_ms2(ms2file):#from experimental
    thres=[]
    with open(ms2file) as ms_sc:
        for nn,line in enumerate(ms_sc):
            if line.rstrip():
                lsp=line.split()
                if lsp[0]=='scan':
                    dp_scan=[]
                else:
                    ms1mz=float(lsp[0])
                    lsp=next(ms_sc).split()
                    rt=float(lsp[0])
                    lsp=next(ms_sc).split()
                    mz_list=[float(x) for x in lsp]
                    lsp=next(ms_sc).split()
                    I_list=[float(x) for x in lsp]

                    I_mz_list=sorted(((y,x) for x,y in zip(mz_list,I_list)),reverse=True)

                    mz_i_list=sorted((y,x) for x,y in I_mz_list[:10])
                    if mz_i_list:
                        i_list=[x for _,x in mz_i_list]
                        dp_scan.append(Spec(ms1mz,rt,[x for x,_ in mz_i_list],i_list,[x for _,x in I_mz_list],[x for x,_ in I_mz_list]))
    return sorted(dp_scan)

def readms1peak(bn):
    ms1peaks=[]
    with open('ms1feature_'+bn+'.txt') as ms1peakfile:
        for line in ms1peakfile:
            lsp=line.rstrip().split()
            if len(lsp)==5:
                ms1peaks.append(Peak(*[float(x) for x in lsp],float(lsp[0])))
    ms1peaks.sort()
    iso_diff=1.00335
    single=0
    doub=0
    peak_double=set()
    for ii in range(len(ms1peaks)):
        peak0=ms1peaks[ii]
        err_bd=min(.01,bound_ppm(peak0.mz*ms1ppm))


        pos0=bisect_left(ms1peaks,(peak0.mz-iso_diff-err_bd,))
        pos1=bisect_left(ms1peaks,(peak0.mz-iso_diff+err_bd,),lo=pos0)
        if pos0!=pos1:
            peak1=min(ms1peaks[pos0:pos1],key=lambda p: abs(peak0.rt-p.rt))
            if abs(peak0.rt-peak1.rt)<1 and peak0.auc<peak1.auc:
                ms1peaks[ii]=Peak(*peak0[:-1],peak1.mmz)
                if peak1.mmz==peak1.mz:
                    single+=1
    return ms1peaks,peak_double

def read_scans(bn):
    ms1scans=[]
    rtset=[]
    with open('ms_scans_'+bn+'.txt') as ms_sc:
        for nn,line in enumerate(ms_sc):
            if line.rstrip():
                lsp=line.rstrip().split(' ',1)
                if lsp[0]=='scan':
                    swath_name=lsp[1]
                    dp_scan=[]
                else:
                    rt=float(lsp[0])
                    rtset.append(rt)
                    lsp=next(ms_sc).rstrip().split()
                    mz_list=[float(x) for x in lsp]
                    lsp=next(ms_sc).rstrip().split()
                    I_list=[float(x) for x in lsp]
                    for mz,I in zip(mz_list,I_list):
                        dp_scan.append(Point(mz,rt,I))
            else:
                if swath_name=='MS1':
                    ms1scans=sorted(dp_scan)
    return ms1scans,sorted(rtset)


def print_score(mzML_file):
    basename0=os.path.basename(mzML_file)
    print(basename0)

    ms2scans=read_ms2('ms2spectra_'+basename0+'.txt')



    ms1peaks,peak_double=readms1peak(basename0)


    ins_f=open('isf_'+basename0+'.txt','w')

    def print_ms2(ms2,p):
        ins_f.write('{} {}\n'.format(p.mz,p.rt))
        ins_f.write(' '.join(str(x) for x in ms2.mz_all)+'\n')
        ins_f.write(' '.join(str(x) for x in ms2.I_all)+'\n')

    def ms2_spectrum():
        ms2_wp=[] #ms2 with ms1 peak info
        for nn,spec in enumerate(ms2scans):
            err_bd=bound_ppm(spec.ms1mz*ms1ppm)
            pos0=bisect_left(ms1peaks,(spec.ms1mz-err_bd,))
            pos1=bisect_left(ms1peaks,(spec.ms1mz+err_bd,),lo=pos0)
            s_peak=None
            if pos0!=pos1:
                peak=min(ms1peaks[pos0:pos1],key=lambda p: abs(spec.rt-p.rt))
                if abs(spec.rt-peak.rt)<peak.sc*1.5:
                    s_peak=peak
            if s_peak is None or s_peak.mz==s_peak.mmz: #append monoisotopic peak only
                ms2_wp.append((spec,s_peak))

        peak_ms2=dict()
        ms2_wp_=[]
        for s,p in ms2_wp:
            if p is None:
                ms2_wp_.append((s,p))
            elif p not in peak_ms2 or abs(s.rt-p.rt)<abs(peak_ms2[p].rt-p.rt):
                peak_ms2[p]=s
        for p,s in peak_ms2.items():
            ms2_wp_.append((s,p))

        ms2_wp_=sorted(ms2_wp_)
        for s_p in ms2_wp_[:]:
            s,_=s_p
            if s_p in ms2_wp_:
                peak_mz=sorted((x for x in ms2_wp_[bisect_left(ms2_wp_,((s.ms1mz,),)):bisect_left(ms2_wp_,((s.ms1mz+.01,),))] if abs(x[0].rt-s.rt)<50),key=lambda x:max(x[0].I),reverse=True)
                for s_p in peak_mz[1:]:
                    ms2_wp_.remove(s_p)

        ms2_wp=sorted(ms2_wp_)

        isf_sc=[[] for x in range(len(ms2_wp))]
        for ii in range(len(ms2_wp)-1,-1,-1):
            ms2,peak=ms2_wp[ii]
            if peak:# and peak.mz==peak.mmz:
                for jj in range(bisect_left(ms2_wp,((ms2.ms1mz,),))):
                    i_ms2,i_p=ms2_wp[jj]
                    if i_p:
                        err_bd=min(.01,bound_ppm(i_p.mz*ms2ppm))
                        pos0=bisect_left(ms2.mz,i_p.mz-err_bd)
                        pos1=bisect_left(ms2.mz,i_p.mz+err_bd,lo=pos0)
                    if i_p and abs(peak.rt-i_p.rt)<ISF_rt_diff and pos0!=pos1 and max(ms2.I)*.1<max(ms2.I[pos0:pos1]): #require ISF to be n% of base peak
                        ms2_I=[]
                        i_ms2_I=[]
                        xfrag=set()
                        for f_i,f_mz in ((x,y) for x,y in zip(i_ms2.I_all[:10],i_ms2.mz_all) if y<i_ms2.ms1mz+.01):
                            err_bd=.01
                            pos0=bisect_left(ms2.mz,f_mz-err_bd)
                            pos1=bisect_left(ms2.mz,f_mz+err_bd,lo=pos0)
                            if pos0!=pos1:
                                i_ms2_I.append(f_i)
                                ms2_I.append(max(ms2.I[pos0:pos1]))
                                for i in range(pos0,pos1): xfrag.add(i)
                            else:
                                i_ms2_I.append(f_i)
                                ms2_I.append(0)
                        if sum(1 for x in ms2_I if x>0)>1:
                            for nn,(f_mz,f_I) in enumerate(zip(ms2.mz,ms2.I)):
                                if nn not in xfrag and f_mz<i_ms2.ms1mz+.01:
                                    ms2_I.append(f_I)
                                    i_ms2_I.append(0)
                            isf_par0=[(x,y) for x,y in zip(i_ms2_I,ms2_I) if y>0]
                            isf_par1=[(x,y) for x,y in zip(i_ms2_I,ms2_I) if x>0]
                            cs=max(cos_sim([x for x,y in isf_par0],[y for x,y in isf_par0]), \
                                    cos_sim([x for x,y in isf_par1],[y for x,y in isf_par1]))
                            if cs>ISF_score:
                                isf_sc[jj].append((cs,ii))

        for jj,ic in enumerate(isf_sc):
            for cs,ii in ic:
                print_ms2(*ms2_wp[ii])
                print_ms2(*ms2_wp[jj])
                ins_f.write('{}\n\n'.format(cs))
        return ms2_wp,isf_sc

    ms2_wp,isf_sc=ms2_spectrum()


    def mass_matching(jj):
        spec,peak=ms2_wp[jj]
        adduct_match=[]

        if "NoMatch" in lib_types[0]:
            for adduct0,(mass0,charge0,_) in list(adduct_list.items())[:-1]:
                Mmass=((spec.ms1mz if peak is None else peak.mz)*charge0-mass0)/(1 if adduct0[0]=='M' else int(adduct0[0]))
                err_bd=bound_ppm(Mmass*ms1ppm)
                pos0=bisect_left(dlist,(Mmass-err_bd,))
                pos1=bisect_left(dlist,(Mmass+err_bd,),lo=pos0)
                adduct_match.append((adduct0,(pos0,pos1)))
        else:
            score_ent=[]
            premz=(spec.ms1mz if peak is None else peak.mz)

            if 1:
                err_bd=bound_ppm(premz*ms1ppm)
                pos_0=bisect_left(lib_ent,(premz-err_bd,))
                pos_1=bisect_left(lib_ent,(premz+err_bd,),lo=pos_0)
                lib_ent_=(x for x in lib_ent[pos_0:pos_1] if x.rt==None or abs(spec.rt-x.rt)<RT_shift)

                for ent in lib_ent_:
                    adduct0,charge0=ent.adduct,ent.charge
                    ms2_I=[]
                    ent_I=[]
                    xfrag=set()
                    for nn,(f_mz,f_I) in enumerate((x,y) for x,y in zip(ent.mz,ent.I) if (charge0*premz-x)>3.3):
                        err_bd=.01
                        pos0=bisect_left(spec.mz,f_mz-err_bd)
                        pos1=bisect_left(spec.mz,f_mz+err_bd,lo=pos0)
                        ent_I.append(f_I)
                        if pos0!=pos1:
                            ms2_I.append(max(spec.I[pos0:pos1]))
                            for i in range(pos0,pos1): xfrag.add(i)
                        else:
                            ms2_I.append(0)
                    m_peaks=sum(1 for x in ms2_I if x>0)
                    if m_peaks>=min_peaks:
                        for nn,(f_mz,f_I) in enumerate(zip(spec.mz,spec.I)):
                            if nn not in xfrag and (charge0*premz-f_mz)>3.3:
                                ms2_I.append(f_I)
                                ent_I.append(0)
                        cs=cos_sim(ent_I,ms2_I)#*sum(x for x,y in zip(ent_I,ms2_I) if y>0)/sum(ent_I)
                        if cs>MS2_score:
                            score_ent.append((adduct0,cs,ent,m_peaks))
            if score_ent:
                max_score_ent=max(score_ent,key=operator.itemgetter(1)) #pick top scoring entry
                adduct_match.append(max_score_ent)

        for cs,ii in isf_sc[jj]:
            _,peak1=ms2_wp[ii]
            adduct_match.append(('?',cs,Ent(0,'ISF of (m/z={:.5f}, rt={:.1f}s) {:.5f}'.format(peak1.mz,peak1.rt,peak.mz),[0],[0],None,None,None,''),0))
        return adduct_match,spec,peak

    ms1scans,rtall=read_scans(basename0)

    def print_ann(ann_,adduct,spec,peak,name):
        ann_.write('NAME:\n')
        ann_.write(name+'\n')
        adduct_c=adduct_list.get(adduct[0],('','',''))
        if adduct[0]=='?': #if isf
            ann_.write('ADDUCT: ISF\n')
        else:
            ann_.write('ADDUCT: {}\n'.format(adduct[0],adduct_c[1],adduct_c[2]))
        ann_.write('TARGET_M/Z, FEATURE_M/Z: {:.6f}'.format(spec.ms1mz))
        ann_.write(', no_ms1_feature_detected' if peak is None else ', '+format(peak.mz,'.6f'))
        ann_.write('\n')

        ann_.write('FORMULA: {}\n'.format(adduct[2].formu))

        ann_.write('SCAN_START_TIME, RT: {:.3f}'.format(spec.rt))

        if peak is None:
            rt_l,rt_r=spec.rt-10,spec.rt+10
            p_area=[x for x in ms1scans[bisect_left(ms1scans,(spec.ms1mz-.01,)):bisect_left(ms1scans,(spec.ms1mz+.01,))] if rt_l<x.rt<rt_r]
        else:
            rt_l,rt_r=peak.rt-peak.sc*1.5,peak.rt+peak.sc*1.5
            p_area=[x for x in ms1scans[bisect_left(ms1scans,(peak.mz-.01,)):bisect_left(ms1scans,(peak.mz+.01,))] if rt_l<x.rt<rt_r]
        ms1rt=rtall[bisect_left(rtall,rt_l):bisect_left(rtall,rt_r)]
        p_dict=dict() # highest intensities per scan bounded by m/z
        for pt in p_area:
            if pt.rt not in p_dict or p_dict[pt.rt]<pt.I:
                p_dict[pt.rt]=pt.I
        p_maxI=[p_dict.get(rt,0.) for rt in ms1rt]
        ms1_auc=sum((I0+I1)*(rt1-rt0) for rt0,rt1,I0,I1 in zip(ms1rt,ms1rt[1:],p_maxI,p_maxI[1:]))/2

        ann_.write(', no_ms1_feature_detected' if peak is None else ', '+format(peak.rt,'.3f'))
        ann_.write('\n')
        if len(adduct)==4:
            dotp=adduct[1]
        elif len(adduct)==2:
            dotp=False
        ann_.write('PEAK_AREA: '+str(ms1_auc)+'\n')
        ann_.write('DOT_PRODUCT: {}\n'.format(format(dotp,'.3f') if dotp else '-'))
        ann_.write('MATCHING_PEAKS: {}\n'.format(adduct[3]))
        ann_.write('EXPERIMENTAL_SPECTRUM:\n')
        for i,mz in zip(spec.I_all,spec.mz_all):
            ann_.write('{:.6f} {:.3f}\n'.format(mz,i))


    def rec_all():
        uk_count=1
        if "NoMatch" not in lib_types[0]:
            una_=open('una_'+basename0+'.txt','w')
        id_quant_ma=collections.defaultdict(list)
        for adduct_match,spec,peak in map(mass_matching, range(len(ms2_wp))):
            for adduct in adduct_match:
                if len(adduct)==4:
                    name0=adduct[2].name
                elif len(adduct)==2:
                    pos0,pos1=adduct[1]
                    nameset={x for _,x in dlist[pos0:pos1]}
                    name0='\n'.join(sorted(nameset))
                id_quant_ma[name0,adduct[0]].append((adduct,spec,peak))
            if not adduct_match and "NoMatch" not in lib_types[0]:
                una_.write("NAME: unknown_{} {} MS1 feature\n".format(uk_count,('with' if peak else 'no')))
                uk_count+=1
                una_.write("PRECURSORMZ: {:.6f}\n".format(spec.ms1mz))
                if ispos:
                    una_.write("PRECURSORTYPE: [unknown]+\n")
                else:
                    una_.write("PRECURSORTYPE: [unknown]-\n")
                una_.write("RETENTIONTIME: {:.3f}\n".format(spec.rt/60))
                thres=bpfilter[0]*max(spec.I_all)
                I_mz_list=[(i,mz) for i,mz in zip(spec.I_all,spec.mz_all) if i>thres]
                I_mz_list=I_mz_list[:bpfilter[1]]
                una_.write("Num Peaks: {}\n".format(len(I_mz_list)))
                for i,mz in I_mz_list:
                    una_.write('{:.6f} {:.3f}\n'.format(mz,i))
                una_.write('\n')
        return id_quant_ma

    id_quant_ma=rec_all()


    def print_all(id_quant):
        annotated_=open('ann_'+basename0+'.txt','w')
        for (name,adductid),adducts in id_quant.items():
            for adduct,spec,peak in adducts:
                print_ann(annotated_,adduct,spec,peak,name)
                if len(adduct)==4:
                    ent=adduct[2]
                    annotated_.write('LIBRARY_SPECTRUM:\n')
                    for mz,i in zip(ent.mz,ent.I):
                        annotated_.write('{} {}\n'.format(mz,i))
                annotated_.write('\n')

    print_all(id_quant_ma)


list(map(print_score, mzML_files))








