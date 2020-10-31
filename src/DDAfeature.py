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
from array import array
import time
import cwt
import commonfn

start_time = time.time()

param_set={
        "mzML_files",
        }

param_dict=commonfn.read_param(param_set)

Point=collections.namedtuple('Point',('rt mz I'))

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


mzML_files=sorted(glob.glob(param_dict["mzML_files"]))



min_group_size=2#int(param_dict["min_group_size"])
mz_space=.009


def print_eic_ms(mzML_file):

    basename0=os.path.basename(mzML_file)
    print(basename0)
    ms1_scans=[]
    ms2_scans=[]

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
                ...
            elif mslevel=='2':
                ms2_scans.append((element.find(".//*[@accession='MS:1000744']").get('value'),store_scan(element)))
            else:
                sys.exit()
        element.clear()
    del element
    del tree

    print(len(ms1_scans),' MS1 scans')
    print(len(ms2_scans),' MS2 scans')


    def mz_slice(ms1_scans,sorted_ms2_pre_mz):
        rtdict={rt:n for n,rt in enumerate(sorted({sc.rt for sc in ms1_scans}))}
        ofile.write('scan MS1\n')
        ofile.write('\t'.join([str(x) for x in rtdict])+'\n')
        data_points=[Point(scan.rt,mz,i) for scan in ms1_scans for mz,i in zip(scan.mz,scan.I)]
        data_points.sort(key=operator.attrgetter('mz'))
        mz_min,mz_max= sorted_ms2_pre_mz[0]-.02,sorted_ms2_pre_mz[-1]+.03

        mzlist=array('d',(mz for _,mz,_ in data_points))
        slice_cut=[]
        for i in itertools.takewhile(lambda n:n<mz_max,itertools.count(mz_min,mz_space)):
            pos = bisect_left(mzlist, i)
            slice_cut.append((pos,i))
        for (pos,mz0),(pos1,mz1) in zip(slice_cut,slice_cut[3:]):
            if pos+min_group_size<pos1:# and max(I for _,_,I in dp_sub)>min_highest_I:
                dp_sub=data_points[pos:pos1]#.tolist()
                ms2pos0=bisect_left(sorted_ms2_pre_mz,mz0+.008999)
                ms2pos1=bisect_left(sorted_ms2_pre_mz,mz1-.008999)
                if ms2pos0<ms2pos1:
                    eic_dict=dict() # highest intensity in this m/z range
                    for rt,mz,I in dp_sub:
                        if rt not in eic_dict or eic_dict[rt][1]<I:
                            eic_dict[rt]=(mz,I)
                    for rt,(mz,i) in sorted(eic_dict.items()):
                        ofile.write('{}\t{}\t{}\n'.format(rt,mz,i))
                    ofile.write('-\n')
        ofile.write('\n')



    with open('eic_'+basename0+'.txt','w') as ofile:
        mz_slice(ms1_scans,sorted(float(x[0]) for x in ms2_scans))


    def print_pt2(ms_scans):
        ofile.write('scan MS2\n')
        for ms1mz,scan_i in ms_scans:
            ofile.write(ms1mz+'\n')
            ofile.write(str(scan_i.rt)+'\n')
            ofile.write(' '.join(str(x) for x in scan_i.mz)+'\n')
            ofile.write(' '.join(str(x) for x in scan_i.I)+'\n')
        ofile.write('\n')

    if len(ms2_scans):
        with open('ms2spectra_'+basename0+'.txt','w') as ofile:
            print_pt2(ms2_scans)

    def print_pt(ms_scans):
        ofile.write('scan MS1\n')
        for scan_i in ms_scans:
            ofile.write(str(scan_i.rt)+'\n')
            ofile.write(' '.join(str(mz) for mz,i in zip(scan_i.mz,scan_i.I) if i>0)+'\n')
            ofile.write(' '.join(str(i) for i in scan_i.I if i>0)+'\n')
        ofile.write('\n')

    with open('ms_scans_'+basename0+'.txt','w') as ofile:
        print_pt(ms1_scans)

    ms2_scans=sorted(ms2_scans,key=lambda x:float(x[0]))
    sorted_ms2_pre_mz=[float(x[0]) for x in ms2_scans]

    basename0=os.path.basename(mzML_file)

    with open('ms1feature_'+basename0+'.txt','w') as writepeak:
        peak_list=[]
        for peaks in map(cwt.findridge, cwt.get_EICs(basename0)):
            peak_list.extend(peaks)
        
        peak_list.sort()
        for peak in peak_list[:]:
            if peak in peak_list:
                peak_mz=sorted((x for x in peak_list[bisect_left(peak_list,(peak.mz,)):bisect_left(peak_list,(peak.mz+.01,))] if abs(x.rt-peak.rt)<x.sc+peak.sc),key=operator.attrgetter('coef'),reverse=True)
                for peak0 in peak_mz[1:]:
                    peak_list.remove(peak0)

        writepeak.write('MS1\n')
        for peak in peak_list:
            pos0=bisect_left(sorted_ms2_pre_mz,peak.mz-.01)
            pos1=bisect_left(sorted_ms2_pre_mz,peak.mz+.01)
            for _,scan_i in ms2_scans[pos0:pos1]:
                if abs(peak.rt-scan_i.rt)<peak.sc:
                    writepeak.write('\t'.join(str(x) for x in peak)+'\n')
                    break


list(map(print_eic_ms, mzML_files))

print("Run time = {:.1f} mins".format(((time.time() - start_time)/60)))


