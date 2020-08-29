import xml.etree.ElementTree as ET
import base64
import struct
import zlib
import sys
import collections
import operator
import itertools
from bisect import bisect_left
import os
import glob
import concurrent.futures
from multiprocessing import freeze_support
from array import array
import time
import cwt
import commonfn

start_time = time.time()

param_set={
        "mzML_files",
        "min_highest_I",
        "num_threads",
        "library",
        "RT_shift",
        "pos/neg mode",
        "ms1_ppm",
        }
param_dict=commonfn.read_param(param_set)

mzML_files=sorted(glob.glob(param_dict["mzML_files"]))

num_threads=int(param_dict["num_threads"])

ispos=True if param_dict['pos/neg mode']=='pos' else False

ms1ppm=float(param_dict['ms1_ppm'])/1e6

min_group_size=2#int(param_dict["min_group_size"])
min_highest_I=float(param_dict["min_highest_I"])
group_I_threshold=min_highest_I#float(param_dict["group_I_threshold"])
lib_types=param_dict["library"].splitlines()
RT_shift=float(param_dict['RT_shift'])

Point=collections.namedtuple('Point',('rt mz I'))
Ent=collections.namedtuple('Ent',('Mmass name adduct rt'))

def bin2float(node):
    d=base64.b64decode(node.findtext("./{http://psi.hupo.org/ms/mzml}binary"))
    if node.find("*/[@accession='MS:1000574']") is not None:
        d=zlib.decompress(d)
    fmt='<{}f'.format(int(len(d)/4)) if node.find("*/[@accession='MS:1000523']") is None else '<{}d'.format(int(len(d)/8))
    return struct.unpack(fmt, d)


def store_scan(element):
    rt=element.find(".//*[@accession='MS:1000016']")
    rtinsec=float(rt.get('value'))
    if rt.get('unitName')=="minute":
        rtinsec*=60
    mz=bin2float(element.find(".//*[@accession='MS:1000514'].."))
    inten=bin2float(element.find(".//*[@accession='MS:1000515'].."))
    return Point(rtinsec,array('d',(m for m,i in zip(mz,inten) if i>0)),array('d',(i for i in inten if i>0)))



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
if 'msdial' in lib_types:
    if ispos:
        libpaths.append(script_dir+'MSMS-Public-Pos-VS12.msp')
    else:
        libpaths.append(script_dir+'MSMS-Public-Neg-VS12.msp')

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


def print_eic_ms(mzML_file):

    basename0=os.path.basename(mzML_file)
    print(basename0)
    ms1_scans=['MS1']
    ms2_scans=['MS2']

    tree=ET.parse(open(mzML_file,'rb'))

    for element in tree.iter(tag='{http://psi.hupo.org/ms/mzml}spectrum'):
        if element.findtext(".//*{http://psi.hupo.org/ms/mzml}binary"):
            mslevel_elem=element.find("*[@accession='MS:1000511']")
            if mslevel_elem is None:
                ref_id=element.find("{http://psi.hupo.org/ms/mzml}referenceableParamGroupRef").attrib['ref']
                mslevel_elem=(tree.find(".//*{http://psi.hupo.org/ms/mzml}referenceableParamGroup[@id='"+ref_id+"']/*[@accession='MS:1000511']"))
                centroid_elem=(tree.find(".//*{http://psi.hupo.org/ms/mzml}referenceableParamGroup[@id='"+ref_id+"']/*[@accession='MS:1000127']"))
            else:
                centroid_elem=element.find("*[@accession='MS:1000127']")
            if centroid_elem is None:
                print("error: profile mode!")
                sys.exit()
            mslevel=mslevel_elem.attrib['value']
            if mslevel=='1':
                ms1_scans.append(store_scan(element))
            elif mslevel=='2':
                ...
            else:
                sys.exit()
        element.clear()
    del element
    del tree

    print(len(ms1_scans)-1,' MS1 scans')


    def mz_slice(ms_scans):
        rtdict={rt:n for n,rt in enumerate(sorted({sc.rt for sc in ms_scans[1:]}))}
        ofile.write('scan '+ms_scans[0]+'\n')
        ofile.write('\t'.join([str(x) for x in rtdict])+'\n')
        data_points=[Point(scan.rt,mz,i) for scan in ms_scans[1:] for mz,i in zip(scan.mz,scan.I)]
        data_points.sort(key=operator.attrgetter('mz'))


        mzlist=array('d',(mz for _,mz,_ in data_points))
        for ent in lib_ent:
            for ii in [-.03,-.015,0]:
                mbd=ent.Mmass+ii
                pos0 = bisect_left(mzlist, mbd)
                pos1 = bisect_left(mzlist, mbd+.03)
                if ent.rt!='NA':
                    dp_sub=[x for x in data_points[pos0:pos1] if abs(ent.rt-x[0])<RT_shift+30]
                else:
                    dp_sub=data_points[pos0:pos1]
                if len(dp_sub)>min_group_size and max(I for _,_,I in dp_sub)>min_highest_I:
                    eic_dict=dict() # highest intensity in this m/z range
                    for rt,mz,I in dp_sub:
                        if rt not in eic_dict or eic_dict[rt][1]<I:
                            eic_dict[rt]=(mz,I)
                    if min_group_size<=len({r for r,(_,i) in eic_dict.items() if i>group_I_threshold}):
                        for rt,mz_i in sorted(eic_dict.items()):
                            ofile.write('{}\t{}\t{}\n'.format(rt,*mz_i))
                        ofile.write('-\n')
        ofile.write('\n')


    with open('eic_'+basename0+'.txt','w') as ofile:
        mz_slice(ms1_scans)

    cwt.cwt(mzML_file)


list(map(print_eic_ms, mzML_files))

import DDArextab
DDArextab.print_tab(lib_ent)

print("Run time = {:.1f} mins".format(((time.time() - start_time)/60)))
