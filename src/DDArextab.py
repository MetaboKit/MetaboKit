import sys
import collections
import operator
import itertools
from bisect import bisect_left
import os
import glob
from array import array
from time import sleep
import bisect
import statistics
import commonfn
from commonfn import bound_ppm


param_set={
        "mzML_files",
        "min_highest_I",
        "library",
        "RT_shift",
        "pos/neg mode",
        "ms1_ppm",
        }
param_dict=commonfn.read_param(param_set)

Peak=collections.namedtuple('Peak',('mz rt sc coef auc'))
Ent=collections.namedtuple('Ent',('Mmass name adduct rt'))

mzML_files=sorted(glob.glob(param_dict["mzML_files"]))

RT_shift=float(param_dict['RT_shift'])
ms1ppm=float(param_dict['ms1_ppm'])/1e6


ms1peaks_dict=dict()
for n,mzML_file in enumerate(mzML_files,1):
    print(n,mzML_file[:-5])
    basename0=os.path.basename(mzML_file)
    ms1peaks=[]
    with open('ms1feature_'+basename0+'.txt') as ms1peakfile:
        for line in ms1peakfile:
            lsp=line.rstrip().split()
            if len(lsp)==5:
                ms1peaks.append(Peak(*[float(x) for x in lsp]))
    ms1peaks_dict[mzML_file]=sorted(ms1peaks)


def print_tab(lib_ent):
    with open('quant_rex.txt','w') as quant_auc:
        quant_auc.write('name\tISF\tadduct\tfeature_m/z(library)\tRT(library)\tfeature_m/z(experimental median)\t%detected\t')
        quant_auc.write('\t'.join(x[:-5] for x in mzML_files)+'\t')
        quant_auc.write('\t'.join('RT_'+x[:-5] for x in mzML_files)+'\n')
        for ent in lib_ent:
            ent_p=[]
            for nn,mzML_file in enumerate(mzML_files):
                ms1peaks=ms1peaks_dict[mzML_file]
                bd=bound_ppm(ent.Mmass*ms1ppm)
                pos0=bisect_left(ms1peaks,(ent.Mmass-bd,))
                pos1=bisect_left(ms1peaks,(ent.Mmass+bd,))
                if ent.rt!='NA':
                    peak=[p for p in ms1peaks[pos0:pos1] if abs(p.rt-ent.rt)<RT_shift]
                else:
                    peak=ms1peaks[pos0:pos1]
                if peak:
                    if ent.rt!='NA':
                        peak=[max(peak,key=operator.attrgetter('auc'))]
                    for p in peak:
                        ent_p.append((nn,p))
            line_str=[]
            if ent_p:
                for nn in range(len(mzML_files)):
                    pos0=bisect_left([x[0] for x in ent_p],nn)
                    pos1=bisect.bisect_right([x[0] for x in ent_p],nn,lo=pos0)
                    line_str.append(( \
                            ','.join(format(dat_n[1].auc,'.1f') for dat_n in ent_p[pos0:pos1]), \
                            ','.join(format(dat_n[1].rt/60,'.2f') for dat_n in ent_p[pos0:pos1]) ))
                mzmed=statistics.median(x[1].mz for x in ent_p)
                if ent.rt!='NA':
                    quant_auc.write('{}\t{}\t{}\t{:.5f}\t{:.2f}\t{:.5f}'.format(ent.name,'*'if ent.name.startswith('ISF of ')else'',ent.adduct,ent.Mmass,ent.rt/60,mzmed))
                else:
                    quant_auc.write('{}\t{}\t{:.5f}\tNA\t{:.5f}'.format(ent.name,ent.adduct,ent.Mmass,mzmed))
                quant_auc.write('\t{:.2f}'.format(sum((1 if x else 0) for x,_ in line_str)/len(mzML_files)))
                quant_auc.write('\t'+'\t'.join(x for x,_ in line_str)+'\t'+'\t'.join(x for _,x in line_str)+'\n')




if __name__ == '__main__':
    ispos=True if param_dict['pos/neg mode']=='pos' else False

    lib_types=param_dict["library"].splitlines()

    libpaths=[]
    script_dir=os.path.abspath(os.path.dirname(sys.argv[0]))+'/libs/'
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
    metabokit=[]
    if any(x.startswith('metabokit ')  for x in lib_types):
        for x in lib_types:
            if x.startswith('metabokit '):
                libpaths.append(x)
    
    def get_cpds():#from libs
        lib_ent=[]
        for libpath in libpaths:
            lib_ent.extend(commonfn.read_msp(libpath))
        lib_ent.sort()
        lib_ent0=[]
        for ent in lib_ent[:]:
            if ent in lib_ent:
                ent_sub=[lib_ent[bisect_left(lib_ent,(ent.Mmass,))]]
                names=[]
                for ent0 in ent_sub:
                    names.append(ent0.name)
                    lib_ent.remove(ent0)
                lib_ent0.append(Ent(ent.Mmass,'---'.join(names),ent.adduct,ent.rt))
        return sorted(lib_ent0)
    
    lib_ent=get_cpds()

    print_tab(lib_ent)
