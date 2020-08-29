import glob
import collections
import os
import itertools
import sys
import re

import statistics
import bisect
import commonfn
from commonfn import bound_ppm


param_set={
        "mzML_files",
        "library",
        "ms1_ppm",
        "ISF_rt_diff",
        "RT_shift",
        "MS2_score",
        "pfcor",
        }
param_dict=commonfn.read_param(param_set)

RT_shift=float(param_dict['RT_shift'])
lib_types=param_dict["library"].splitlines()
ms1ppm=float(param_dict['ms1_ppm'])/1e6
ISF_rt_diff=float(param_dict["ISF_rt_diff"])
MS2_score=float(param_dict['MS2_score'])
pfcor=float(param_dict['pfcor'])


count_id=collections.Counter()
mzML_files=sorted(glob.glob(param_dict["mzML_files"]))
quantfiles=['quant_'+'_'.join(lib_types)+'_'+os.path.basename(x)+'.txt' for x in mzML_files]
for n,mzML_file in enumerate(mzML_files,1):
    print(n,mzML_file[:-5])


merge_q=collections.defaultdict(lambda:[('','','','')]*len(mzML_files))
merge_q_=collections.defaultdict(lambda:['']*len(mzML_files))#precursor quant before isf adding
merge_isf=collections.defaultdict(lambda:[('','','','')]*len(mzML_files))

for quantfile in quantfiles:
    with open(quantfile) as qf:
        qf.readline()
        for lsp in (line.rstrip().split('\t') for line in qf):
            if lsp[1]=='MS1' and float(lsp[6])>0 and float(lsp[7])>-1:
                count_id[lsp[0]+lsp[2]]+=1

for ss,quantfile in enumerate(quantfiles):
    with open(quantfile) as qf:
        line1=qf.readline()
        for lsp in (line.rstrip().split('\t') for line in qf):
            if count_id[lsp[0]+lsp[2]]>0.0*len(mzML_files):
                if False:
                    merge_isf[tuple(float(x) for x in re.findall(r'\d+\.\d+',lsp[0]))+(0 if lsp[1]=='MS1' else float(lsp[3]),)][ss]=(lsp[-1],lsp[5],lsp[6],lsp[7])
                else:
                    merge_q[tuple(lsp[:5])][ss]=[lsp[-1],lsp[5],lsp[6],lsp[7]]
                    if lsp[1]=='MS1' and lsp[-1]:
                        merge_q_[tuple(lsp[:3])][ss]=float(lsp[-1])

merge_frag=collections.defaultdict(list)
merge_prec=collections.defaultdict(list)
for cpd,dat in merge_q.items():
    if cpd[1]!='MS1':
        merge_frag[cpd[0],cpd[2]].append((cpd,[(float(x) if x else 0) for x,_,_,_ in dat]))
    else:
        merge_prec[cpd[0],cpd[2]]=(cpd,[(float(x) if x else 0) for x,_,_,_ in dat])


sort_isf=sorted(merge_isf.keys())
merge_qp=[(x,y) for x,y in merge_q.items() if x[1]=='MS1']
for cpd,dat in merge_qp:
    if not re.fullmatch("\d+\.\d+",cpd[4]): break
    pmz,rt=[float(x) for x in cpd[3:5]]
    pos0=bisect.bisect_left(sort_isf,(pmz-pmz*ms1ppm,))
    pos1=bisect.bisect_left(sort_isf,(pmz+pmz*ms1ppm,))
    for isf in sort_isf[pos0:pos1]:
        if abs(rt-isf[1])<ISF_rt_diff:
            if isf[3]==0:#if precursor
                for nn,((isf_q,isf_rt,_,_),(par_q,par_rt,_,_)) in enumerate(zip(merge_isf[isf],dat)):
                    if isf_q and par_q:
                        if abs(float(isf_rt)-float(par_rt))>3:#ISF_rt_diff:
                            print(nn,cpd,isf_rt,par_rt)
                        merge_q[cpd][nn][0]=str(float(par_q)+float(isf_q))

merge_cpd=collections.defaultdict(list)
for cpd,dat in merge_q.items():
    merge_cpd[cpd[0]+cpd[2]].append((cpd,dat))

p_mz_rt=[]
for name_,dat in (x for x in merge_q.items() if x[0][1]=='MS1'):
    if name_[0].startswith('ISF of'):
        pmz=float(name_[0].split('m/z=')[1].split(',')[0])
        rt=float(name_[0].split('rt=' )[1].split(')')[0])
    else:
        pmz=float(name_[3])
        rt=(float(name_[4]) if re.fullmatch("\d+\.\d+",name_[4]) else None)
    pos0=bisect.bisect_left(p_mz_rt,(pmz-.01,))
    pos1=bisect.bisect_left(p_mz_rt,(pmz+.01,),lo=pos0)
    for x in p_mz_rt[pos0:pos1]:
        if rt is None or abs(rt-x[1])<RT_shift:
            break
    else:
        bisect.insort(p_mz_rt,(pmz,rt))
adduct_list=set()
mz_rt_name=collections.defaultdict(list)
for name_,dat in (x for x in merge_q.items() if x[0][1]=='MS1'):
    if name_[0].startswith('ISF of'):
        pmz=float(name_[0].split('m/z=')[1].split(',')[0])
        rt=float(name_[0].split('rt=' )[1].split(')')[0])
    else:
        pmz=float(name_[3])
        rt=(float(name_[4]) if re.fullmatch("\d+\.\d+",name_[4]) else None)
    pos0=bisect.bisect_left(p_mz_rt,(pmz-.01,))
    pos1=bisect.bisect_left(p_mz_rt,(pmz+.01,),lo=pos0)
    for n,x in enumerate(p_mz_rt[pos0:pos1],pos0):
        if rt is None or abs(rt-x[1])<RT_shift:
            mz_rt_name[n].append(name_[0]+'\n'+name_[2])#name+adduct
            adduct_list.add(name_[2])
            break
    else:
        print('err')
        sys.exit()

grp_count=1
name_mz_rt=dict()
name_isf_grp=dict()
for key,names in sorted(mz_rt_name.items()):
    if any((not x.startswith('ISF of')) for x in names):
        for name in names:
            name_mz_rt[name]=grp_count
        grp_count+=1
        isf_pres=any(x.startswith('ISF of') for x in names)#check for ISF presence
        for name in names:
            name_isf_grp[name]=isf_pres
mz_rt_n_dict={k:collections.defaultdict(list) for k in list(adduct_list)}#+['isf']}
mz_rt_n=[]
for k,n in name_mz_rt.items():
    if k.startswith('ISF of'):
        mz_rt_n.append((n,[k]))
    else:
        mz_rt_n_dict[k.rsplit('\n',1)[1].split()[0]][n].append(k)

pf_dict=collections.defaultdict(list)
for cpd,dat in merge_q.items():
    pf_dict[cpd[0]+'\n'+cpd[2]].append(cpd)

n0_names=dict() #rename ISFs
for n0,names in sorted([inner for outer in [list(x.items()) for x in mz_rt_n_dict.values()] for inner in outer]+mz_rt_n):
    if not names[0].startswith('ISF of '):
        n0_names[n0]=[x.split('\n')[0] for x in names]

def format_fn(vec,fmt_str):
    return [format(float(x),fmt_str) if x else '' for x in vec]

with open('quant_'+'_'.join(lib_types)+'All.txt','w') as mq:
    lsp=line1.rstrip().split('\t')
    headx=[x[:-5] for x in mzML_files]
    mq.write('group\t'+'\t'.join(lsp[:5])+'\t'+'\t'.join(headx)+'\t'+'\t'.join('RT_'+x for x in headx)+'\t'+'\t'.join('DOTP_'+x for x in headx)+'\t'+'\t'.join('COR_'+x for x in headx)+'\n')


    for n0,names in sorted([inner for outer in [list(x.items()) for x in mz_rt_n_dict.values()] for inner in outer]+mz_rt_n):
        for name in names:

            for cpd in pf_dict[name]:# for precursor and frags in the precursor grp
                if cpd[1]=='MS1':
                    dat=merge_q[cpd]
                    countdot=sum(float(x)>MS2_score  for _,_,x,_ in dat)
                    countcor=sum(float(x)>pfcor  for _,_,_,x in dat)
            if countdot<3 or countcor<3:
                continue

            for cpd in pf_dict[name]:# for precursor and frags in the precursor grp
                dat=merge_q[cpd]
                quant_=[x for x,_,_,_ in dat]
                if cpd[1]=='MS1':
                    pre_quant=[float(x) if x else x for x,_,_,_ in dat]
                else:
                    for nn,(par_isf,par,f) in enumerate(zip(pre_quant,merge_q_[(cpd[0],'MS1',cpd[2])],quant_)):
                        if par and par_isf and f and par_isf>par:#isf adjustment for frag 
                            quant_[nn]=str(float(f)*par_isf/par)#frag_quant*((pre_quant+isf_quant)/pre_quant)
                cpd_name='ISF of '+'---'.join(n0_names[n0]) if cpd[0].startswith('ISF of ') else cpd[0]
                mq.write(str(n0)+'\t'+cpd_name+'\t'+'\t'.join(cpd[1:])+'\t'+'\t'.join(format_fn(quant_,'.2f'))+'\t'+'\t'.join(format_fn([x[1] for x in dat],'.1f'))+'\t'+'\t'.join(format_fn([x[2] for x in dat],'.2f'))+'\t'+'\t'.join(format_fn([x[3] for x in dat],'.2f'))+'\n')
                

