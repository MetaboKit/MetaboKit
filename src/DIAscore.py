
import collections
import operator
import sys
import math
import statistics
from bisect import bisect_left
import re
import os
import pathlib
import glob
import concurrent.futures
import xml.etree.ElementTree as ET
from multiprocessing import freeze_support

import commonfn
from commonfn import bound_ppm


param_set={
        "mzML_files",
        "library",
        "ms1_ppm",
        "ms2_ppm",
        "num_threads",
        "rt_diff",
        "topNfrag",
        "MS2_score",
        "pfcor",
        "ms2_auc_w/o_feature",
        }
param_dict=commonfn.read_param(param_set)

Point=collections.namedtuple('Point',('mz rt I'))
Peak=collections.namedtuple('Peak',('mz rt sc coef auc mmz'))

lib_types=param_dict["library"].splitlines()

ms1ppm=float(param_dict['ms1_ppm'])/1e6
ms2ppm=float(param_dict['ms2_ppm'])/1e6

num_threads=int(param_dict["num_threads"])
rt_diff=int(param_dict["rt_diff"])
topNfrag=int(param_dict["topNfrag"])

MS2_score=float(param_dict['MS2_score'])
pfcor=float(param_dict['pfcor'])

mzML_files=sorted(glob.glob(param_dict["mzML_files"]))

ms2_auc_no_feat=int(param_dict["ms2_auc_w/o_feature"])

libpaths=[]
script_dir=os.path.abspath(os.path.dirname(sys.argv[0]))+'/libs/'
if 'nist' in lib_types:
    libpaths.append(script_dir+'nist.Original.bak.msp')
if 'LipidBlast' in lib_types:
    libpaths.append(script_dir+'LipidBlast-ASCII-spectra/LipidBlast-neg.msp')
    libpaths.extend(glob.glob(script_dir+'LipidBlast-ASCII-spectra/custom-libs/*neg.msp'))
    libpaths.append(script_dir+'LipidBlast-ASCII-spectra/LipidBlast-pos.msp')
    libpaths.extend(glob.glob(script_dir+'LipidBlast-ASCII-spectra/custom-libs/*pos.msp'))
if 'LipidBlast-fork' in lib_types:
    libpaths.append(script_dir+'MSDIAL-InsilicoMSMS-Lipids-Neg.msp')
    libpaths.append(script_dir+'MSDIAL-InsilicoMSMS-Lipids-Pos.msp')
if 'hmdb' in lib_types:
    libpaths.append(script_dir+'hmdb_metabolites.xml')
if 'msdial' in lib_types:
    libpaths.append(script_dir+'MSMS-Public-Pos-VS11.msp')
    libpaths.append(script_dir+'MSMS-Public-Neg-VS11.msp')
if 'NoMatch' in lib_types[0]:
    libpaths=[script_dir+'Database_Dec2017.txt']
    if len(lib_types[0])>9:
        libpaths=[lib_types[0].split(' ',1)[1]]
        print(lib_types,libpaths)
if 'metabokit' in lib_types:
    for file0 in glob.glob('ann_*All.txt'):
        libpaths.append(file0)
print(libpaths)

cpd_list={
        'CO' :27.99491, 'H'  :1.007825, 'Li' :7.016004, 'NH4':18.03437, 'Na' :22.98977, 'Na2':45.97954,
        'H'  :1.007825, '2H' :2.015650,
        '2H' :2.015650, '2I' :253.8089, '2K' :77.92741, '2Na':45.97954, '3H' :3.023475, '3K' :116.8911, '3Na':68.96931, 'H'  :1.007825, 'H2O':18.01056, 'I'  :126.9045, 'K'  :38.96371, 'NH3':17.02655, 'NH4':18.03437, 'Na' :22.98977, 'OH' :17.00274, 
        }

def cos_sim(list1,list2):
    if len(list1)!=len(list2):
        print('adf')
        sys.exit()
    if sum(list1)<=0 or sum(list2)<=0: return 0
    return sum(x*y for x,y in zip(list1,list2))/math.sqrt(sum(x*x for x in list1)*sum(x*x for x in list2))

def cos_sim0(list1,list2):
    if len(list1)!=len(list2):
        print('adf')
        sys.exit()
    list1=list1[:topNfrag]
    list2=list2[:topNfrag]
    if sum(list1)<=0 or sum(list2)<=0: return 0
    return sum(math.sqrt(x*y) for x,y in zip(list1,list2))/math.sqrt(sum(list1)*sum(list2))
    return sum(x*y for x,y in zip(list1,list2))/math.sqrt(sum(x*x for x in list1)*sum(x*x for x in list2))

def read_lib(libpath):
    adductset=set()
    lib_dict=collections.defaultdict(list)

    if "LipidBlast" in libpath and "All.txt" not in libpath:
        it_cpd=iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            lsp=line.split(': ')
            if lsp[0]=="Name":
                name=lsp[1]
                adduct=line[line.find("[M")+1:line.rfind("]")]
                adduct=tuple(re.split('(\-|\+)',adduct))
                charge=line[line.rfind("]")+1:]
                charge=charge[1] if charge[0]=='(' else charge[0]
                charge=(charge if charge.isdigit() else '1')
            elif lsp[0]=="PRECURSORMZ":
                ms1mz=(lsp[1])
            elif lsp[0]=="Comment":
                name+=' '+lsp[1].rstrip(';').rsplit('; ',1)[-1]
            elif lsp[0]=="Num Peaks":
                frag_mz=[]
                frag_I=[]
                for lsp in (line.split() for line in it_cpd):
                    if not lsp: break
                    frag_mz.append((lsp[0]))
                    frag_I.append((lsp[1]))
                lib_dict[ms1mz+' '+charge+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' NA'].append(name)
    if "nist" in libpath and "All.txt" not in libpath:
        it_cpd=iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            lsp=line.split(': ')
            if lsp[0]=="Name":
                name=lsp[1]
                line=next(it_cpd)
                name+=' '+line[line.find("["):]
                adduct=line[line.find("[")+1:line.rfind("]")]
                adduct=tuple(re.split('(\-|\+)',adduct))
                charge=line[line.rfind("]")+1]
                charge=(charge if charge.isdigit() else '1')
            elif lsp[0]=="Formula":
                name+=' '+lsp[1]
            elif lsp[0]=="NISTNO":
                name+=' NIST '+lsp[1]
            elif lsp[0]=="PrecursorMZ":
                ms1mz=(lsp[1])
            elif lsp[0]=="Num peaks":
                frag_mz=[]
                frag_I=[]
                for lsp in (line.split() for line in it_cpd):
                    if not lsp: break
                    frag_mz.append((lsp[0]))
                    frag_I.append((lsp[1]))
                if re.fullmatch("\d+\.\d{2,}",frag_mz[0]) is None: continue #remove if <2dp
                if charge!=0 and all((x in cpd_list) for x in adduct[2::2]):
                    lib_dict[ms1mz+' '+charge+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' NA'].append(name)
    if 'hmdb' in libpath and "All.txt" not in libpath:
        hmdb_dict=dict()
        for _, elem in ET.iterparse(open(libpath,'rb')):
            if elem.tag=='{http://www.hmdb.ca}metabolite':
                accession=elem.findtext('{http://www.hmdb.ca}accession')
                name=elem.findtext('{http://www.hmdb.ca}name')
                chem_for=elem.findtext('{http://www.hmdb.ca}chemical_formula')
                mono_weight=elem.findtext('{http://www.hmdb.ca}monisotopic_molecular_weight')
                if re.fullmatch("\d+\.\d{2,}",mono_weight):
                    hmdb_dict[accession]=(name,chem_for,mono_weight)
                elem.clear()
        for file_x in (glob.glob(str(pathlib.Path.home())+'/AMD_lib/hmdb_experimental_msms_spectra/*')):
            frag_mz=[]
            frag_I=[]
            tree=ET.parse(open(file_x,'rb'))
            hmdb_id=tree.findtext('database-id')
            for peak in tree.iter(tag='ms-ms-peak'):
                frag_mz.append(peak.find('mass-charge').text)
                frag_I.append(peak.find('intensity').text)
            hmdb_dat=hmdb_dict.get(hmdb_id)
            datbase=tree.findtext('./references/reference/database')
            datbase=('unknown' if datbase is None else datbase)
            if hmdb_dat is not None and frag_mz:
                lib_dict[str(float(hmdb_dat[2])+cpd_list['H'])+' 1 '+','.join(frag_mz)+' '+','.join(frag_I)+' '+'M,+,H'+' NA'].append(hmdb_dat[0]+' '+hmdb_dat[1]+' '+datbase+' '+tree.findtext('id'))
    if 'VS12' in libpath and "All.txt" not in libpath:
        it_cpd=iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            if line.startswith('NAME: '):
                name=line.split(': ')[1]
            elif line.startswith('PRECURSORMZ: '):
                ms1mz=line.split(': ')[1]
            elif line.startswith('PRECURSORTYPE: '):
                adduct=line.split(': ')[1][:-1].strip('[]')
                adductset.add(adduct)
                adduct=tuple(re.split('(\-|\+)',adduct))
            elif line.startswith('FORMULA: '):
                name+=' '+line.split(': ')[1]
            elif line.startswith('Num Peaks: '):
                frag_mz=[]
                frag_I=[]
                for lsp in (line.split() for line in it_cpd):
                    if not lsp: break
                    frag_mz.append((lsp[0]))
                    frag_I.append((lsp[1]))
                if frag_mz and all((x in cpd_list) for x in adduct[2::2]):
                    lib_dict[ms1mz+' 1 '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' NA'].append(name)
    if "All.txt" in libpath:
        it_ann=iter(line.rstrip() for line in open(libpath))
        for line in it_ann:
            if line.startswith("NAME:"):
                name=next(it_ann)
                for line in it_ann:
                    if not line.startswith("ADDUCT: "):
                        name+='---'+line
                    else:
                        adduct_charge=line.split(': ')[1].split()
                        if len(adduct_charge)==1:
                            adduct=adduct_charge
                            charge='1+'
                        else:
                            *adduct,charge=adduct_charge
                        break
                for line in it_ann:
                    if line.startswith("PRECURSOR_M/Z: "):
                        ms1mz=line.split(': ')[1]
                        break
                line=next(it_ann) #rt
                rt=line.split(': ')[1]
                for line in it_ann:
                    if line.startswith("LIBRARY_SPECTRUM:"):
                        break
                frag_mz=[]
                frag_I=[]
                frag_ann=[]
                for line in it_ann:
                    if line:
                        mz,I=line.split(' ')[:2]
                        frag_mz.append(mz)
                        frag_I.append(I)
                        if len(line.split(' '))>2:
                            frag_ann.append(line.split(' ',2)[2])
                        else:
                            frag_ann.append('MS2')
                    else:
                        break
                lib_dict[ms1mz+' '+charge[:-1]+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' '+rt+' '+','.join(frag_ann)].append(name)
    return lib_dict

def get_cpds():
    lib_dict=dict()
    for libpath in libpaths:
        lib_dict.update(read_lib(libpath))
    return lib_dict.items()

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
    return ms1peaks,peak_double

def read_scans(bn):
    ms1scans=[]
    ms2scans=dict()
    with open('ms_scans_'+bn+'.txt') as ms_sc:
        for nn,line in enumerate(ms_sc):
            if line.rstrip():
                lsp=line.rstrip().split(' ',1)
                if lsp[0]=='scan':
                    swath_name=lsp[1]
                    dp_scan=[]
                else:
                    rt=float(lsp[0])
                    lsp=next(ms_sc).rstrip().split()
                    mz_list=[float(x) for x in lsp]
                    lsp=next(ms_sc).rstrip().split()
                    I_list=[float(x) for x in lsp]
                    for mz,I in zip(mz_list,I_list):
                        dp_scan.append(Point(mz,rt,I))
            else:
                if swath_name=='MS1':
                    ms1scans=sorted(dp_scan)
                else:
                    ms2scans[swath_name]=sorted(dp_scan)

    sswath=sorted(ms2scans.keys(),key=lambda x: float(x.split()[0]))
    startpt=[float(x.split()[0]) for x in sswath[1:]]
    end__pt=[float(x.split()[1]) for x in sswath[:-1]]
    minmz=float(sswath[0].split()[0])
    maxmz=float(sswath[-1].split()[1])
    return ms1scans,ms2scans,sswath,startpt,end__pt,minmz,maxmz


def print_score(mzML_file):
    basename0=os.path.basename(mzML_file)
    print(basename0)

    ms1peaks,peak_double=readms1peak(basename0)


    for n,line in enumerate(open('eic_'+basename0+'.txt')):
        if n==1:
            lsp=line.rstrip().split()
            rtdict={float(rt):n for n,rt in enumerate(lsp)}
            break
    rtset=sorted(rtdict.keys())


    ms1scans,ms2scans,sswath,startpt,end__pt,minmz,maxmz=read_scans(basename0)

    cpd_quant=open('quant_'+'_'.join(lib_types)+'_'+basename0+'.txt','w')
    cpd_quant.write('compound\tQuantification Mode\tadduct\tmass\tRT(library)\trt\tdot_prod\tp_f_cor\tquant\n')


    def scoring(cpd):
        ent,name=cpd
        ms1mz,charge,frag_mz,frag_I,adduct,RT,frag_ann=ent.split(' ',6)
        ms1mz=float(ms1mz)
        if re.fullmatch("\d+\.\d+",RT): RT=float(RT)
        frag_I=[float(x) for x in frag_I.split(',')]
        frag_mz=[float(x) for x in frag_mz.split(',')]
        frag_ann=[x for x in frag_ann.split(',')]
        sorted_I=sorted(zip(frag_I,frag_mz),reverse=True)[:topNfrag]
        frag_mz=[x for _,x in sorted_I]
        frag_I=[x for x,_ in sorted_I]
        if len(frag_I)!=len(frag_mz): print('abb')
        adduct=adduct.replace(',','')
        name='---'.join(name)
        
        frag_mz_l=[0]*len(frag_mz)
        frag_mz_r=[0]*len(frag_mz)
        for nn,f_mz in enumerate(frag_mz):
            err_bd=.01
            frag_mz_l[nn]=f_mz-err_bd
            frag_mz_r[nn]=f_mz+err_bd

        err_bd=bound_ppm(ms1mz*ms1ppm)
        pos0=bisect_left(ms1peaks,(ms1mz-err_bd,))
        pos1=bisect_left(ms1peaks,(ms1mz+err_bd,))
        pseudo_feat=[]#pseudo feature for entries w/o feature
        ms1peaks_match=[x for x in ms1peaks[pos0:pos1] if not isinstance(RT,float) or abs(RT-x.rt)<rt_diff]
        if ms2_auc_no_feat and not ms1peaks_match and isinstance(RT,float):#if no feature found and RT in library
            pseudo_feat.append(Peak(mz=ms1mz,rt=RT,sc=10,coef=0,auc=0,mmz=ms1mz))
        score_peaks=[]
        for ms1peak in ms1peaks_match+pseudo_feat:
            if not minmz<ms1peak.mz<maxmz:
                continue
            pos0=bisect_left(startpt,ms1peak.mz)
            pos1=bisect_left(end__pt,ms1peak.mz)
            if pos0==pos1: # ms1peak in one window only
                pos=pos0
            else: # take the window whose boundaries are furthest from ms1peak
                pos=(pos0 if ms1peak.mz-startpt[pos0-1]>end__pt[pos1]-ms1peak.mz else pos1)
            iso=ms2scans[sswath[pos]]
            ms2_I=[0]*len(frag_I)
            ms2_auc=[0]*len(frag_I)
            pfc=[0]*len(frag_I)
            rt_l=ms1peak.rt-ms1peak.sc*1.5
            rt_r=ms1peak.rt+ms1peak.sc*1.5
            ms1rt=rtset[bisect_left(rtset,rt_l):bisect_left(rtset,rt_r)]
            p_dict=dict() # highest intensities per scan bounded by m/z
            p_area=[x for x in ms1scans[bisect_left(ms1scans,(ms1peak.mz-.01,)):bisect_left(ms1scans,(ms1peak.mz+.01,))] if rt_l<x.rt<rt_r]
            for pt in p_area:
                if pt.rt not in p_dict or p_dict[pt.rt]<pt.I:
                    p_dict[pt.rt]=pt.I
            p_maxI=[p_dict.get(rt,0.) for rt in ms1rt]
            for nn,(f_mz_l,f_mz_r,f_I) in enumerate(zip(frag_mz_l,frag_mz_r,frag_I)):
                f_area=[x for x in iso[bisect_left(iso,(f_mz_l,)):bisect_left(iso,(f_mz_r,))] if rt_l<x.rt<rt_r]
                if f_area:
                    f_area_=[x.I for x in f_area if abs(x.rt-ms1peak.rt)<2]
                    if f_area_:
                        ms2_I[nn]=max(f_area_)
                        f_dict=dict() # highest intensities per scan bounded by m/z
                        for pt in f_area:
                            if pt.rt not in f_dict or f_dict[pt.rt]<pt.I:
                                f_dict[pt.rt]=pt.I
                        ms2_maxI=[f_dict.get(rt,0.) for rt in ms1rt]
                        pfc[nn]=cos_sim(p_maxI,ms2_maxI)
                        ms2_auc[nn]=sum((I0+I1)*(rt1-rt0) for rt0,rt1,I0,I1 in zip(ms1rt,ms1rt[1:],ms2_maxI,ms2_maxI[1:]))/2
            ms1_auc=sum((I0+I1)*(rt1-rt0) for rt0,rt1,I0,I1 in zip(ms1rt,ms1rt[1:],p_maxI,p_maxI[1:]))/2
            ssm=cos_sim0(frag_I,ms2_I)
            score_peaks.append((ssm,pfc,ms1peak,ms2_auc,ms1_auc))
        if score_peaks:
            score_peaks.sort(reverse=True)
            max_score_peaks=[score_peaks[0]]
            for x in score_peaks[1:]:
                if max_score_peaks[0][0]-x[0]<.1:
                    max_score_peaks.append(x)
                else:
                    break
            max_peak=max(max_score_peaks,key=lambda x:x[2].auc) #pick top scoring ms2 for each entry, if score difference is insignificant use auc
            return (name,adduct,ms1mz,RT,frag_mz,frag_ann,max_peak)
        return False


    def print_mp(mp):
        name,adduct,ms1mz,RT,frag_mz,frag_ann,max_peak=mp
        if len(max_peak[1])>2:
            ms1pfc=sorted(max_peak[1])[-1]#nth largest pfc
        else:
            ms1pfc=min(max_peak[1])
        if True:
            cpd_quant.write('{}\tMS1\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name,adduct,ms1mz,RT,max_peak[2].rt,max_peak[0],ms1pfc,max_peak[4]))
            for f_mz,f_ann,ms2_pfc,ms2_auc in zip(frag_mz,frag_ann,max_peak[1],max_peak[3]):
                if ms2_pfc>-.5:
                    cpd_quant.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(name,f_ann,adduct,f_mz,RT,max_peak[2].rt,max_peak[0],ms2_pfc,ms2_auc))


    for mp in map(scoring,get_cpds()):
        if mp: print_mp(mp)


list(map(print_score, mzML_files))



