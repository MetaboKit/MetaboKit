
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
import commonfn
from commonfn import bound_ppm

param_set={
        "mzML_files",
        "library",
        "adduct",
        "ms1_ppm",
        "ms2_ppm",
        "background",
        "ISF_rt_diff",
        "ISF_score",
        "MS2_score",
        "min_peaks",
        "num_threads",
        "base_peak_filter",
        "RT_shift",
        }
param_dict=commonfn.read_param(param_set)

Point=collections.namedtuple('Point',('mz rt I'))
Spec=collections.namedtuple('Spec',('ms1mz rt mz I mz_all I_all'))
Peak=collections.namedtuple('Peak',('mz rt sc coef auc mmz'))
Ent=collections.namedtuple('Ent',('Mmass name mz I adduct charge rt'))

RT_shift=float(param_dict['RT_shift'])
lib_types=param_dict["library"].splitlines()

ms1ppm=float(param_dict['ms1_ppm'])/1e6
ms2ppm=float(param_dict['ms2_ppm'])/1e6

num_threads=int(param_dict["num_threads"])
min_peaks=int(param_dict["min_peaks"])
bpfilter=param_dict.get("base_peak_filter",'0.00 20').split(' ')
bpfilter=[float(bpfilter[0]),int(bpfilter[1])]

ISF_rt_diff=float(param_dict["ISF_rt_diff"])
ISF_score=float(param_dict["ISF_score"])
MS2_score=float(param_dict["MS2_score"])

adduct_list={x.split()[0][0:]:(float(x.split()[1]),int(x.split()[2][0]),x.split()[2][1]) for x in param_dict["adduct"].splitlines()}
pos_neg_set={x[-1] for x in adduct_list.values()}
if pos_neg_set=={'+'}:
    ispos=True
elif pos_neg_set=={'-'}:
    ispos=False
else:
    print(list(adduct_list.items()))
    print('check adduct list')
    sys.exit()

adduct_list['?']=(0, 1, '?')
print(list(adduct_list.items())[:-1])

mzML_files=sorted(glob.glob(param_dict["mzML_files"]))
libpaths=[]
script_dir=os.path.abspath(os.path.dirname(sys.argv[0]))+'/libs/'
if 'nist' in lib_types:
    libpaths.append(script_dir+'nist.Original.bak.msp')
if 'LipidBlast' in lib_types:
    if ispos:
        libpaths.append(script_dir+'LipidBlast-ASCII-spectra/LipidBlast-pos.msp')
        libpaths.extend(glob.glob(script_dir+'LipidBlast-ASCII-spectra/custom-libs/*pos.msp'))
    else:
        libpaths.append(script_dir+'LipidBlast-ASCII-spectra/LipidBlast-neg.msp')
        libpaths.extend(glob.glob(script_dir+'LipidBlast-ASCII-spectra/custom-libs/*neg.msp'))
if 'LipidBlast-fork' in lib_types:
    if ispos:
        libpaths.append(script_dir+'MSDIAL-InsilicoMSMS-Lipids-Pos.msp')
    else:
        libpaths.append(script_dir+'MSDIAL-InsilicoMSMS-Lipids-Neg.msp')
if 'hmdb' in lib_types:
    libpaths.append(script_dir+'hmdb_metabolites.xml')
if 'msdial' in lib_types:
    if ispos:
        libpaths.append(script_dir+'MSMS-Public-Pos-VS12.msp')
    else:
        libpaths.append(script_dir+'MSMS-Public-Neg-VS12.msp')
if 'NoMatch' in lib_types[0]:
    libpaths=[script_dir+'Database_Dec2017.txt']
    if len(lib_types[0])>9:
        libpaths=[lib_types[0].split(' ',1)[1]]
        print(lib_types,libpaths)

metabokit=[]
if any(x.startswith('metabokit ')  for x in lib_types):
    for x in lib_types:
        if x.startswith('metabokit '):
            libpaths.append(x)



cpd_list={
        'CO' :27.99491, 'H'  :1.007825, 'Li' :7.016004, 'NH4':18.03437, 'Na' :22.98977, 'Na2':45.97954,
        'H'  :1.007825, '2H' :2.015650,
        '2H' :2.015650, '2I' :253.8089, '2K' :77.92741, '2Na':45.97954, '3H' :3.023475, '3K' :116.8911, '3Na':68.96931, 'H'  :1.007825, 'H2O':18.01056, 'I'  :126.9045, 'K'  :38.96371, 'NH3':17.02655, 'NH4':18.03437, 'Na' :22.98977, 'OH' :17.00274, 
        '':0,
        'FA':46.00548, 'Hac':60.02113,
        'ACN':41.02655,
        }

def cos_sim(list1,list2):
    if len(list1)!=len(list2):
        print('adf')
        sys.exit()
    if sum(list1)<=0 or sum(list2)<=0: return 0
    return sum(math.sqrt(x*y) for x,y in zip(list1,list2))/math.sqrt(sum(list1)*sum(list2))
    return sum(x*y for x,y in zip(list1,list2))/math.sqrt(sum(x*x for x in list1)*sum(x*x for x in list2))


def read_lib(libpath):
    adductset=set()
    lib_dict=collections.defaultdict(list)
    libpath0=libpath.split('/libs/')[-1]
    if 'LipidCreatorValidStudy_MRM_Workklist_V1_forHyungWon.txt' in libpath0:
        qq=open(libpath)
        qq.readline()
        nnnn=0
        for nnn,line in enumerate(qq):
            lsp=line.split('\t')
            name="Bo's_MRM_list "+lsp[1]+' '+lsp[11]
            ms1mz=lsp[2]
            adduct=lsp[11]
            charge=adduct[adduct.rfind(']')+1]
            charge=(charge if charge.isdigit() else '1')
            adduct=adduct[adduct.find("[")+1:adduct.rfind("]")].replace(' ','')
            adductset.add(adduct)
            adduct=tuple(x.strip() for x in re.split('(\-|\+)',adduct))

            frag_mz=[lsp[3]]
            frag_I=['1']
            if frag_mz:
                lib_dict[ms1mz+' '+charge+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' NA'].append(name)
                nnnn+=1
        print(list(lib_dict)[0])
        print(adductset)
        print(nnn,nnnn)
    if "LipidBlast-ASCII-spectra" in libpath0:
        it_cpd=iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            lsp=line.split(': ',1)
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
                    frag_mz.append(lsp[0])
                    frag_I.append(lsp[1])
                lib_dict[ms1mz+' '+charge+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' NA'].append(name)

    if 'MSDIAL-InsilicoMSMS-Lipids-' in libpath0:
        it_cpd=iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            if line.startswith('NAME: '):
                name=line.split(': ',1)[1]
            elif line.startswith('PRECURSORMZ: '):
                ms1mz=line.split(': ',1)[1]
            elif line.startswith('PRECURSORTYPE: '):
                adduct=line.split(': ',1)[1]
                adduct=adduct[adduct.find('[')+1:adduct.rfind(']')]
                adduct=tuple(re.split('(\-|\+)',adduct))
                charge=line[line.rfind("]")+1:]
                charge=charge[1] if charge[0]=='(' else charge[0]
                charge=(charge if charge.isdigit() else '1')
            elif line.startswith('FORMULA: '):
                name+=' '+line.split(': ',1)[1]
            elif line.startswith('Num Peaks: '):
                frag_mz=[]
                frag_I=[]
                for lsp in (line.split() for line in it_cpd):
                    if not lsp: break
                    frag_mz.append(lsp[0])
                    frag_I.append(lsp[1])
                lib_dict[ms1mz+' '+charge+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' NA'].append(name)

    if "nist.Original.bak.msp" in libpath0:
        it_cpd=iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            lsp=line.split(': ',1)
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
                name+=' NIST '#+lsp[1]
            elif lsp[0]=="PrecursorMZ":
                ms1mz=(lsp[1])
            elif lsp[0]=="Num peaks":
                frag_mz=[]
                frag_I=[]
                for lsp in (line.split() for line in it_cpd):
                    if not lsp: break
                    frag_mz.append(lsp[0])
                    frag_I.append(lsp[1])
                if re.fullmatch("\d+\.\d{2,}",frag_mz[0]) is None: continue #remove if <n dp
                lib_dict[ms1mz+' '+charge+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' NA'].append(name)
    if '-VS12.msp' in libpath0 or libpath0.startswith('metabokit '):
        if libpath0.startswith('metabokit '):
            libpath=libpath0.split(' ',1)[1]
            if not open(libpath):
                print(libpath+" not found")
                sys.exit()
        it_cpd=iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            if line.startswith('NAME: '):
                name=line.split(': ',1)[1]
                RT='NA'
            elif line.startswith('PRECURSORMZ: '):
                ms1mz=line.split(': ',1)[1]
            elif line.startswith('PRECURSORTYPE: '):
                adduct=line[line.find("[")+1:line.rfind("]")]
                adduct=tuple(re.split('(\-|\+)',adduct))
                charge=line[line.rfind("]")+1]
                charge=(charge if charge.isdigit() else '1')
            elif line.startswith('FORMULA: '):
                name+=' '+line.split(': ',1)[1]
            elif line.startswith('RETENTIONTIME: ') and libpath0.startswith('metabokit '):
                RT=str(float(line.split(': ',1)[1])*60)
            elif line.startswith('Num Peaks: '):
                frag_mz=[]
                frag_I=[]
                for lsp in (line.split() for line in it_cpd):
                    if not lsp: break
                    frag_mz.append(lsp[0])
                    frag_I.append(lsp[1])
                if frag_mz:
                    lib_dict[ms1mz+' '+charge+' '+','.join(frag_mz)+' '+','.join(frag_I)+' '+','.join(adduct)+' '+RT].append(name)
    if 'hmdb_metabolites.xml' in libpath0:
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
        for file_x in glob.glob(script_dir+'hmdb_experimental_msms_spectra/*'):
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
                lib_dict[hmdb_dat[2]+' 1 '+','.join(frag_mz)+' '+','.join(frag_I)+' '+'M,+'+' NA'].append(hmdb_dat[0]+' '+hmdb_dat[1]+' '+datbase+' '+tree.findtext('id'))
    print(libpath0)
    print(len(lib_dict))
    return lib_dict


def get_cpds():#from libs
    lib_dict=dict()
    for libpath in libpaths:
        lib_dict.update(read_lib(libpath))
    lib_ent=[]
    for ent,name in lib_dict.items():
        ms1mz,charge,frag_mz,frag_I,adduct,RT=ent.split(' ')
        ms1mz,charge=float(ms1mz),int(charge)
        RT=(None if RT=='NA' else float(RT))
        frag_mz=[float(x) for x in frag_mz.split(',')]
        frag_I=[float(x) for x in frag_I.split(',')]

        I_mz_list=sorted(zip(frag_I,frag_mz),reverse=True)[:10]
        frag_mz=[x for _,x in I_mz_list]
        frag_I=[x for x,_ in I_mz_list]

        name='\n'.join(name)

        lib_ent.append(Ent(ms1mz,name,frag_mz,frag_I,adduct.replace(',',''),charge,RT))
        continue

        Mmass=ms1mz*charge
        adduct=adduct.split(',')
        for pm,cpd in zip(adduct[1::2],adduct[2::2]):
            if pm=='+':
                Mmass-=cpd_list[cpd]
            elif pm=='-':
                Mmass+=cpd_list[cpd]
            else:
                sys.exit()
        if adduct[0][0].isdigit(): #M,2M,3M
            Mmass/=int(adduct[0][0])


        lib_ent.append(Ent(Mmass,name,frag_mz,frag_I,None,None,RT))
    return sorted(lib_ent)


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
    lib_ent=get_cpds()

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

    bg=param_dict.get("background")
    if bg is not None:
        ms2scans_bg=read_ms2('ms2spectra_'+os.path.basename(glob.glob(bg)[0])+'.txt')
        for ii,ms2sc in enumerate(ms2scans):
            err_bd=bound_ppm(ms2sc.ms1mz*ms1ppm)
            pos0=bisect_left(ms2scans_bg,(ms2sc.ms1mz-err_bd,))
            pos1=bisect_left(ms2scans_bg,(ms2sc.ms1mz+err_bd,),lo=pos0)
            if pos1!=pos0: 
                closest_bg=min(ms2scans_bg[pos0:pos1],key=lambda x:abs(x.rt-ms2sc.rt))
                if abs(closest_bg.rt-ms2sc.rt)<60:
                    new_mz=[]
                    new_I=[]
                    for m,i in zip(ms2sc.mz,ms2sc.I):
                        pos0=bisect_left(closest_bg.mz,m-err_bd)
                        pos1=bisect_left(closest_bg.mz,m+err_bd,lo=pos0)
                        if pos0==pos1:
                            new_mz.append(m)
                            new_I.append(i)
                    ms2scans[ii]=Spec(*ms2sc[:2],new_mz,new_I,*ms2sc[-2:])


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

            for iiiiii in range(1):
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
                    if sum(1 for x in ms2_I if x>0)>=min_peaks: #min number of matching peak
                        for nn,(f_mz,f_I) in enumerate(zip(spec.mz,spec.I)):
                            if nn not in xfrag and (charge0*premz-f_mz)>3.3:
                                ms2_I.append(f_I)
                                ent_I.append(0)
                        cs=cos_sim(ent_I,ms2_I)#*sum(x for x,y in zip(ent_I,ms2_I) if y>0)/sum(ent_I)
                        if cs>MS2_score:
                            score_ent.append((adduct0,cs,ent))
            if score_ent:
                max_score_ent=max(score_ent,key=operator.itemgetter(1)) #pick top scoring entry
                adduct_match.append(max_score_ent)

        for cs,ii in isf_sc[jj]: #if isf, attach data
            _,peak1=ms2_wp[ii]
            adduct_match.append(('?',cs,Ent(0,'ISF of (m/z={:.6f}, rt={:.3f}) {:.6f}'.format(peak1.mz,peak1.rt,peak.mz),[0],[0],None,None,None)))
        return adduct_match,spec,peak

    ms1scans,rtall=read_scans(basename0)

    def print_ann(ann_,adduct,spec,peak,name):
        ann_.write('NAME:\n')
        ann_.write(name+'\n')
        adduct_c=adduct_list.get(adduct[0],('','',''))
        if adduct[0]=='?': #if isf
            ann_.write('ADDUCT: -\n')
        else:
            ann_.write('ADDUCT: {} {}{}\n'.format(adduct[0],adduct_c[1],adduct_c[2]))
        ann_.write('TARGET_M/Z, FEATURE_M/Z: {:.6f}'.format(spec.ms1mz))
        ann_.write(', no_ms1_feature_detected' if peak is None else ', '+format(peak.mz,'.6f'))
        ann_.write('\n')
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
        if len(adduct)==3:
            dotp=adduct[1]
        elif len(adduct)==2:
            dotp=False
        ann_.write('PEAK_AREA: '+str(ms1_auc)+'\n')
        ann_.write('DOT_PRODUCT: '+(format(dotp,'.3f') if dotp else '-')+'\n')
        ann_.write('EXPERIMENTAL_SPECTRUM:\n')
        for i,mz in zip(spec.I_all,spec.mz_all):
            ann_.write('{:.6f} {:.6g}\n'.format(mz,i))


    def rec_all():
        uk_count=1
        if "NoMatch" not in lib_types[0]:
            una_=open('una_'+'_'.join(lib_types)+'_'+basename0+'.txt','w')
        id_quant_ma=collections.defaultdict(list)
        for adduct_match,spec,peak in map(mass_matching, range(len(ms2_wp))):
            for adduct in adduct_match:
                if len(adduct)==3: # if lipidblast
                    name0=adduct[2].name
                elif len(adduct)==2: #if nomatch
                    pos0,pos1=adduct[1]
                    nameset={x for _,x in dlist[pos0:pos1]}
                    name0='\n'.join(sorted(nameset))
                id_quant_ma[name0,adduct[0]].append((adduct,spec,peak))
            if not adduct_match and "NoMatch" not in lib_types[0]:
                una_.write("NAME: unknown_{} {} MS1 feature\n".format(uk_count,('with' if peak else 'no')))
                uk_count+=1
                una_.write("PRECURSORMZ: {:.6f}\n".format(spec.ms1mz))
                if ispos:
                    una_.write("PRECURSORTYPE: [M+H]+ [M+2H]2+ [M+Na]+ [M+NH4]+ [M-H2O+H]+\n")
                else:
                    una_.write("PRECURSORTYPE: [M-H]- [M-2H]2- [M-H2O-H]- [M+HCOO]- [M+Cl]-\n")
                una_.write("RETENTIONTIME: {:.3f}\n".format(spec.rt/60))
                thres=bpfilter[0]*max(spec.I_all)
                I_mz_list=[(i,mz) for i,mz in zip(spec.I_all,spec.mz_all) if i>thres]
                I_mz_list=I_mz_list[:bpfilter[1]]
                una_.write("Num Peaks: {}\n".format(len(I_mz_list)))
                for i,mz in I_mz_list:
                    una_.write('{:.6f} {:.6g}\n'.format(mz,i))
                una_.write('\n')
        return id_quant_ma

    id_quant_ma=rec_all()


    def print_all(id_quant):
        annotated_=open('ann_'+'_'.join(lib_types)+'_'+basename0+'.txt','w')
        for (name,adductid),adducts in id_quant.items():
            for adduct,spec,peak in adducts:
                print_ann(annotated_,adduct,spec,peak,name)
                if len(adduct)==3:
                    ent=adduct[2]
                    annotated_.write('LIBRARY_SPECTRUM:\n')
                    for mz,i in zip(ent.mz,ent.I):
                        annotated_.write('{} {}\n'.format(mz,i))
                annotated_.write('\n')

    print_all(id_quant_ma)


list(map(print_score, mzML_files))






