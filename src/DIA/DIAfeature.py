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
        "window_setting",
        "min_highest_I",
        "num_threads",
        }
param_dict=commonfn.read_param(param_set)

Point=collections.namedtuple('Point',('rt mz I'))

def bin2float(node):
    if node is None: return ()
    binary_txt=node.findtext("./{http://psi.hupo.org/ms/mzml}binary")
    if binary_txt=='': return ()
    d=base64.b64decode(binary_txt)
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


mzML_files=sorted(glob.glob(param_dict["mzML_files"]))

num_threads=int(param_dict["num_threads"])


SWATHs=param_dict["window_setting"].splitlines()

min_group_size=2#int(param_dict["min_group_size"])
min_highest_I=float(param_dict["min_highest_I"])
group_I_threshold=min_highest_I#float(param_dict["group_I_threshold"])
mz_space=.015

def print_eic_ms(mzML_file):

    basename0=os.path.basename(mzML_file)
    print(basename0)

    ms1_scans=['MS1']
    ms2_scans=[[x] for x in SWATHs]

    tree=ET.parse(open(mzML_file,'rb'))

    swath_i=len(SWATHs)
    for element in tree.iter(tag='{http://psi.hupo.org/ms/mzml}spectrum'):
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
            if swath_i!=len(SWATHs):
                print('spectrum index',element.get('index'),len(SWATHs),swath_i)
                for ii in range(min(swath_i,len(SWATHs))): ms2_scans[ii].pop()
            swath_i=0
            ms1_scans.append(store_scan(element))
        elif mslevel=='2':
            if swath_i<len(SWATHs):
                ms2_scans[swath_i].append(Point(ms1_scans[-1].rt,*store_scan(element)[1:])) #use ms1 rt
            swath_i+=1
        else:
            sys.exit()
        element.clear()
    del element
    del tree

    print(len(ms1_scans)-1,' MS1 scans')





    def mz_slice(ms_scans):
        rtdict={rt:n for n,rt in enumerate(sorted({sc.rt for sc in ms_scans[1:]}))}
        ofile.write('scan '+ms_scans[0]+'\n')
        ofile.write('\t'.join([str(x) for x in rtdict.keys()])+'\n')
        data_points=[Point(scan.rt,mz,i) for scan in ms_scans[1:] for mz,i in zip(scan.mz,scan.I)]
        data_points.sort(key=operator.attrgetter('mz'))
        mz_min,mz_max=data_points[0].mz,data_points[-1].mz

        mzlist=array('d',(mz for _,mz,_ in data_points))
        slice_cut=[]
        for i in itertools.takewhile(lambda n:n<mz_max,itertools.count(mz_min,mz_space)):
            pos = bisect_left(mzlist, i)
            slice_cut.append(pos)
        slice_cut.append(len(data_points))
        for pos,pos1 in zip(slice_cut,slice_cut[2:]):
            dp_sub=data_points[pos:pos1]#.tolist()
            if pos+min_group_size<pos1 and max(I for _,_,I in dp_sub)>min_highest_I:
                eic_dict=dict() # highest intensity in this m/z range
                for rt,mz,I in dp_sub:
                    if rt not in eic_dict or eic_dict[rt][1]<I:
                        eic_dict[rt]=(mz,I)
                if min_group_size<=len({r for r,(_,i) in eic_dict.items() if i>group_I_threshold}):
                    for rt,(mz,i) in sorted(eic_dict.items()):
                        ofile.write('{}\t{}\t{}\n'.format(rt,mz,i))
                    ofile.write('-\n')
        ofile.write('\n')

                    

    with open('eic_'+basename0+'.txt','w') as ofile:
        mz_slice(ms1_scans)


    def print_pt(ms_scans):
        ofile.write('scan '+ms_scans[0]+'\n')
        for scan_i in ms_scans[1:]:
            ofile.write(str(scan_i.rt)+'\n')
            ofile.write(' '.join(str(mz) for mz,i in zip(scan_i.mz,scan_i.I) if i>0)+'\n')
            ofile.write(' '.join(str(i) for i in scan_i.I if i>0)+'\n')
        ofile.write('\n')

    with open('ms_scans_'+basename0+'.txt','w') as ofile:
        list(map(print_pt, [ms1_scans]+ms2_scans))


list(map(print_eic_ms, mzML_files))
if __name__ == '__main__':
    freeze_support()
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        list(executor.map(cwt.cwt, mzML_files))


print("Run time = {:.1f} mins".format(((time.time() - start_time)/60)))


