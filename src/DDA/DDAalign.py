import glob
import collections
import statistics
import operator
import re
import sys
import os
from bisect import bisect_left
import bisect
import commonfn
from commonfn import bound_ppm


param_set={
        "mzML_files",
        "library",
        "ms1_ppm",
        "ms2_ppm",
        "adduct",
        "RT_shift",
        }
param_dict=commonfn.read_param(param_set)
#count_id=collections.Counter()
RT_shift=float(param_dict['RT_shift'])
lib_types=param_dict["library"].splitlines()
ms1ppm=float(param_dict['ms1_ppm'])/1e6
adduct_list={x.split()[0][0:]:(float(x.split()[1]),int(x.split()[2][0]),x.split()[2][1]) for x in param_dict["adduct"].splitlines()}
mzML_files=sorted(glob.glob(param_dict["mzML_files"]))
ann_files=['ann_'+'_'.join(lib_types)+'_'+os.path.basename(x)+'.txt' for x in mzML_files]
for n,mzML_file in enumerate(mzML_files,1):
    print(n,mzML_file[:-5])

####count occurance of cpd
#for ann_file in ann_files:
#    it_ann=iter(line.rstrip() for line in open(ann_file))
#    name_ids=set()
#    for line in it_ann:
#        if line.startswith("NAME:"):
#            name_id=next(it_ann)
#            for line in it_ann:
#                if not line.startswith("ADDUCT: "):
#                    name_id+='\n'+line #add all putative
#                else:
#                    name_id+='\n'+line.split(': ')[1] # add adduct goto next
#                    break
#            next(it_ann) #precursor_mz
#            next(it_ann) #rt
#            #if float(next(it_ann).split(': ')[1])<.5: #want at least 1 with dotp>.5
#            #    continue
#            if name_id not in name_ids:#count once per sample
#                count_id[name_id]+=1
#                name_ids.add(name_id)


all_dat=collections.defaultdict(list)
lib_dict=dict()
#all_id=set()
for nn,ann_file in enumerate(ann_files):
    it_ann=iter(line.rstrip() for line in open(ann_file))
    for line in it_ann:
        if line.startswith("NAME:"):
            name_id=next(it_ann)
            for line in it_ann:
                if not line.startswith("ADDUCT: "):
                    name_id+='\n'+line
                else:
                    name_id+='\n'+line.split(': ')[1]
                    break
            #if count_id[name_id]>0: #.5*len(ann_files):
            line=next(it_ann)
            #premz=float(line[line.find(': ')+2:line.rfind(', ')])
            premz=line[line.find(': ')+2:].split(', ')
            premz=float(premz[1]) if re.fullmatch("\d+\.\d+",premz[1]) else float(premz[0])
            line=next(it_ann)
            rt=line[line.find(': ')+2:].split(', ')
            rt=float(rt[1]) if re.fullmatch("\d+\.\d+",rt[1]) else float(rt[0])
            line=next(it_ann)#peak_area
            auc=line.split(': ')[1]
            auc=float(auc) if re.fullmatch("\d+\.\d+",auc) else None
            line=next(it_ann)#dotp
            dotp=line.split(': ')[1]
            dotp=float(dotp) if re.fullmatch("\d+\.\d+",dotp) else None
            next(it_ann)#experimental_spec
            mz_list=[]
            I_list=[]
            for line in it_ann:
                if line and not line.startswith("LIBRARY_SPECTRUM:"):
                    mz,I=line.split()
                    mz_list.append(float(mz))
                    I_list.append(float(I))
                else:
                    break
            lib_dat=''
            for line in it_ann:
                if line:
                    lib_dat+=line+'\n'
                else:
                    break
            mz_list=[float(x) for x in mz_list]
            I_list=[x/max(I_list)*999. for x in I_list]
            #all_dat[name_id,lib_dat].append((nn,dotp,premz,rt,mz_list,I_list,auc))
            all_dat[name_id].append((nn,dotp,premz,rt,mz_list,I_list,auc))
            lib_dict[name_id]=(dotp,lib_dat)
            #all_id.add(name_id)

##merge by rt and m/z
#print(list(all_dat.items())[1][0])
#for i in list(all_dat.items())[1][1]:
#    print(i)
#####

#merge cpd with similar -
#merge isf with same parent
isf_dict=dict()
#name_dict=collections.defaultdict(list)
for name_id in all_dat.keys():
    name_add=name_id.split('\n')
    name,adduct=name_add[:-1],name_add[-1]
    name0=name[0]
    if name0.startswith('ISF of'):
        pmz,rt,mz=[float(x) for x in re.findall("\d+\.\d+",name0)[:3]]
        for pmz0,rt0,mz0,adduct0 in isf_dict.keys():
            if abs(pmz-pmz0)<bound_ppm(pmz*ms1ppm) and abs(rt-rt0)<RT_shift and abs(mz-mz0)<bound_ppm(mz*ms1ppm) and adduct==adduct0:
                isf_dict[(pmz0,rt0,mz0,adduct)].append(name_id)
                break
        else:
            isf_dict[(pmz,rt,mz,adduct)]=[name_id]



name_dict0=dict()
#join ISF that were merged
for key,vv in isf_dict.items():
    if len(vv)>1:
        pmz_rt_mz=[tuple(float(x) for x in re.findall("\d+\.\d+",v)[:3]) for v in vv]
        key0=statistics.median(x for x,_,_ in pmz_rt_mz)
        key1=statistics.median(x for _,x,_ in pmz_rt_mz)
        key2=statistics.median(x for _,_,x in pmz_rt_mz)
        isf_name='ISF of (m/z={}, rt={}) {}\n{}'.format(key0,key1,key2,key[3])
        for v in vv:
            name_dict0[v]=isf_name

#name_dict=collections.defaultdict(set)
#for name_id,dat in all_dat.items():
#    name_add=name_id.split('\n')
#    name,adduct=name_add[:-1],name_add[-1]
#    name0=name[0]
#    if not name0.startswith('ISF of'):
#        rt=statistics.median(x[3] for x in dat)
#        mz=statistics.median(x[2] for x in dat)
#        for isf_name in name_dict0.values():
#            rt0,mz0=[float(x) for x in re.findall("\d+\.\d+",isf_name)[1:3]]
#            if abs(rt-rt0)<10 and abs(mz-mz0)<mz*ms1ppm:
#                name_dict[isf_name].add(name_id.split('\n')[0])


all_dat0=collections.defaultdict(list)
lib_dict0=collections.defaultdict(list)
for name,dat in all_dat.items():
    #if name in name_dict0:
    #    name_id=name.split('\n')[0]+'\n'+'\n'.join(list(name_dict[name]))+'\n'+name.split('\n')[1]
    #    print(name)
        
    all_dat0[name_dict0.get(name,name)].extend(dat)
    lib_dict0[name_dict0.get(name,name)].append(lib_dict[name])

#grouping compounds or ISF
p_mz_rt=[]
for name,dat in all_dat0.items():
    if name.startswith('ISF of'):
        pmz,rt=[float(x) for x in re.findall("\d+\.\d+",name)[:2]]
    else:
        pmz=dat[0][2]
        rt=statistics.median(rt for _,_,_,rt,_,_,_ in dat)
    pos0=bisect_left(p_mz_rt,(pmz-.01,))
    pos1=bisect_left(p_mz_rt,(pmz+.01,),lo=pos0)
    #if "C18H30O3S" in name:
    #    print(name,[x[:4] for x in dat])
    #    #sys.exit()
    for x in p_mz_rt[pos0:pos1]:
        if abs(rt-x[1])<RT_shift:
            break
    else:
        bisect.insort(p_mz_rt,(pmz,rt))
        #p_mz_rt.append((pmz,rt))

mz_rt_name=collections.defaultdict(list)
#for name,dat in sorted(all_dat0.items(),key=lambda x:x[1][0][2]):
for name,dat in all_dat0.items():
    if name.startswith('ISF of'):
        pmz,rt=[float(x) for x in re.findall("\d+\.\d+",name)[:2]]
    else:
        pmz=dat[0][2]
        rt=statistics.median(rt for _,_,_,rt,_,_,_ in dat)
    pos0=bisect_left(p_mz_rt,(pmz-.01,))
    pos1=bisect_left(p_mz_rt,(pmz+.01,),lo=pos0)
    for n,x in enumerate(p_mz_rt[pos0:pos1],pos0):
        if abs(rt-x[1])<RT_shift:
            #name_mz_rt[name]=n
            mz_rt_name[n].append(name)
            break
    else:
        print('err')
        sys.exit()


isf_mz_rt=[]
for key,names in sorted(mz_rt_name.items()):
    if any((not x.startswith('ISF of')) for x in names):
        for name in names:
            if name.startswith('ISF of'):
                pmz,rt,mz=[float(x) for x in re.findall("\d+\.\d+",name)[:3]]
                isf_mz_rt.append((mz,rt,pmz))
isf_mz_rt.sort()

#isf_mz_rt=[tuple(float(x) for x in re.findall("\d+\.\d+",isf_name)[:3]) for x in mz_rt_name.values() for isf_name in x if isf_name.startswith('ISF of')]

grp_count=1
name_mz_rt=dict()
name_isf_grp=dict()
possible_isf=set()
#only include ISF with identified parents
for key,names in sorted(mz_rt_name.items()):
    if any((not x.startswith('ISF of')) for x in names):
        
        ###remove compounds with ISF that are possibly ISF
        flag=False
        for name in (x for x in names if not x.startswith('ISF of')):
            if name not in all_dat0: print('daf')
            mz0=statistics.median(x[2] for x in all_dat0[name])
            rt0=statistics.median(x[3] for x in all_dat0[name])
            pos0=bisect_left(isf_mz_rt,(mz0-bound_ppm(mz0*ms1ppm),))
            pos1=bisect_left(isf_mz_rt,(mz0+bound_ppm(mz0*ms1ppm),))
            for mz,rt,_ in isf_mz_rt[pos0:pos1]:
                if abs(rt-rt0)<RT_shift:
                    flag=True
                    possible_isf.add(name)
                    break
        isf_in_names=[x for x in names if x.startswith('ISF of')]
        if flag and len(isf_in_names)>-1:#remove "possible isf" with more than n ISF
            #print(names)
            continue
        #######################

        for name in names:
            name_mz_rt[name]=grp_count
        grp_count+=1
        isf_pres=any(x.startswith('ISF of') for x in names)#check for ISF presence
        #mark groups with ISF
        for name in names:
            name_isf_grp[name]=isf_pres

#print('possible_isf',possible_isf)


#flip over key value
#dict of dict
mz_rt_n_dict={k:collections.defaultdict(list) for k in list(adduct_list.keys())}#+['isf']}
#separate isf and compound
mz_rt_n=[]
for k,n in name_mz_rt.items():
    if k.startswith('ISF of'):
        mz_rt_n.append((n,[k]))
    else:
        ##############using original adduct
        #mz_rt_n_dict[list(mz_rt_n_dict.keys())[0]][n].append(k)
        ##############

        ##############using all adduct
        mz_rt_n_dict[k.rsplit('\n',1)[1].split()[0]][n].append(k)
        ##############

n0_names=dict() #rename ISFs
for n0,names in sorted([inner for outer in [list(x.items()) for x in mz_rt_n_dict.values()] for inner in outer]+mz_rt_n):
    #print(n0,names)
    if not names[0].startswith('ISF of '):
        n0_names[n0]=[x.split('\n')[0] for x in names]
        
#for k,v in (n0_names.items()):
#    print(v)
##mark peak that are ISF
isf_dat=set()
for n0,names in sorted([inner for outer in [list(x.items()) for x in mz_rt_n_dict.values()] for inner in outer]+mz_rt_n):
    if names[0].startswith('ISF of '):
        for nn,dotp,premz,rt,mz_l,I_l,auc in all_dat0[names[0]]:
            isf_dat.add((nn,premz,rt))
#################
#for idat in isf_dat:
#    if 137.03<idat[2]<137.05:
#        print(idat)


with open('ann_'+'_'.join(lib_types)+'All.txt','w') as cpd_ann, \
        open('quant_'+'_'.join(lib_types)+'All.txt','w') as quant_auc:
    quant_auc.write('group\tISF\tname\tadduct\tfeature_m/z\t')
    quant_auc.write('\t'.join(x[:-5] for x in mzML_files)+'\t')
    quant_auc.write('\t'.join('RT_'+x[:-5] for x in mzML_files)+'\t')
    quant_auc.write('\t'.join('score_'+x[:-5] for x in mzML_files)+'\n')
    n1=1
    prev_n0=1
    for n0,names in sorted([inner for outer in [list(x.items()) for x in mz_rt_n_dict.values()] for inner in outer]+mz_rt_n):
        dat=[]
        dat0=collections.defaultdict(list)
        lib_dict1=[]
        for name in names:
            auc_dfdict=collections.defaultdict(list)
            ##remove peaks identified as ISF
            for nn,dotp,premz,rt,mz_l,I_l,auc in all_dat0[name]:
                if name.startswith('ISF of ') or (nn,premz,rt) not in isf_dat:
                    auc_dfdict[nn].append((dotp,premz,rt,mz_l,I_l,auc))
                    lib_dict1.extend(lib_dict0[name])
            ###################
            for nn in range(len(mzML_files)):
                dotp_rt_auc=auc_dfdict.get(nn)

                if dotp_rt_auc:
                    dotp_rt_auc.sort(reverse=True)
                    max_score_ents=[x[:-1]+(x[-1] if x[-1] else 0,) for x in dotp_rt_auc]
                    if len([x for x in max_score_ents if x[-1]])>0:
                        #all top scoring ms2 with ms1 feature
                        dotp_rt_auc=[x for x in max_score_ents if x[-1]]
                    else:#if multiple peak all no feature, select top scoring
                        dotp_rt_auc=max_score_ents[:1]
                    for x in dotp_rt_auc:
                        dat0[nn].append(x)
        if not dat0:
            continue
        for nn,dotp_l in dat0.items():
            sdotp=sorted(dotp_l,reverse=True)
            dat.append([nn]+list(sdotp[0]))
            for p in sdotp[1:]:
                if sdotp[0][0]-p[0]<.1:
                    dat.append([nn]+list(p))
                else:
                    break
        dat.sort()
        #if 195.051<dat[0][2]<195.053:
        #    for iii in dat:
        #        print(iii)
        nameset=sorted(set(x.rsplit('\n',1)[0] for name in names for x in name.split('\n')[:-1]))
        if any((name in possible_isf) for name in names):
            nameset=['possibly an ISF']+nameset
        adductset=sorted(set(name.rsplit('\n',1)[1] for name in names))

        premz_l=[]
        rt_l=[]
        for _,_,premz,rt,_,_,_ in dat:
            premz_l.append(premz)
            rt_l.append(rt)
        cpd_ann.write("NAME:\n")
        cpd_ann.write('\n'.join(nameset)+'\n')
        cpd_ann.write("ADDUCT: "+', '.join(adductset)+'\n')
        cpd_ann.write("SAMPLE, RT, DOT_PRODUCT, PEAK_AREA\n")
        for nn,dotp,premz,rt,mz_l,I_l,auc in dat:
            #cpd_ann.write(mzML_files[nn][-55:-5].ljust(52)+format(rt,'.2f').ljust(9)+format(dotp,'.2f').ljust(6)+(format(auc,'.1f') if auc else 'no_ms1_feature')+"\n")
            cpd_ann.write('{:52}{:<9.2f}{:<6.2f}{:.1f}\n'.format(mzML_files[nn][-55:-5],rt,dotp,auc if auc else 'no_ms1_feature'))
        cpd_ann.write("PRECURSOR_M/Z: {:.5f}\n".format(statistics.median(premz_l)))
        cpd_ann.write("RT: {:.2f}\n".format(statistics.median(rt_l)))
        cpd_ann.write("EXPERIMENTAL_SPECTRUM:\n")
        max_dat=max(dat,key=operator.itemgetter(1))#print the highest scoring spectrum
        for  mz,I in zip(max_dat[4],max_dat[5]):#max_dat[4:6]:
            cpd_ann.write('{:.5f} {:.5g}\n'.format(mz,I))
        cpd_ann.write("LIBRARY_SPECTRUM:\n")
        if len(lib_dict1)==1:
            lib_dat=lib_dict1[0][1]
            cpd_ann.write(lib_dat)
        else:
            ###take entry with the best dotp
            #lib_dat=min(lib_dict1,key=lambda x:x.count('\n'))[1]
            lib_dat=max(lib_dict1)[1]
            cpd_ann.write(lib_dat)
        cpd_ann.write('\n')

        #quant_auc.write(str(name_mz_rt[name])+'\t')
        if prev_n0!=n0: n1+=1
        quant_auc.write(str(n1)+'\t')
        quant_auc.write(('*' if name_isf_grp[name] else '')+'\t')
        if nameset[0].startswith('ISF of '):
            quant_auc.write('ISF of '+'---'.join(n0_names[n0]))
        else:
            quant_auc.write('---'.join(nameset))
        prev_n0=n0
        quant_auc.write('\t{}\t{:.5f}'.format(','.join(adductset),statistics.median(premz_l)))

        #quant_auc.write('\t'+'\t'.join(format(auc,'.1f') if auc else 'no_ms1_feature' for nn,dotp,premz,rt,mz_l,I_l,auc in dat))
        #    quant_auc.write('\t')
        line_str=[]
        for nn in range(len(mzML_files)):
            pos0=bisect_left([x[0] for x in dat],nn)
            pos1=bisect.bisect_right([x[0] for x in dat],nn,lo=pos0)
            #pos1=bisect_left([x[0] for x in dat],nn+1,lo=pos0)
            line_str.append(( \
                    ','.join(format(dat_n[-1],'.1f') if dat_n[-1] else 'no_ms1_feature' for dat_n in dat[pos0:pos1]), \
                    ','.join(format(dat_n[3]/60,'.2f') for dat_n in dat[pos0:pos1]), \
                    ','.join(format(dat_n[1],'.2f') for dat_n in dat[pos0:pos1]) ))
        quant_auc.write('\t'+'\t'.join(x for x,_,_ in line_str)+'\t'+'\t'.join(x for _,x,_ in line_str)+'\t'+'\t'.join(x for _,_,x in line_str)+'\n')


print('done')