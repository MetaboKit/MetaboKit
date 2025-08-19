import glob
import collections
import statistics
import operator
import re
import sys
import os
from bisect import bisect_left
import bisect

import DDAreadlib
import DDAcommonfn
from DDAcommonfn import bound_ppm


param_set = {
    "mzML_files",
    "library",
    "ms1_ppm",
    "ms2_ppm",
    "RT_shift",
}
param_dict = DDAcommonfn.read_param(param_set)
RT_shift = float(param_dict["RT_shift"])
lib_types = [x.split(" #", 1)[0].strip() for x in param_dict["library"].splitlines()]
ms1ppm = float(param_dict["ms1_ppm"]) / 1e6
adduct_list = {"": ""}
mzML_files = sorted(glob.glob(param_dict["mzML_files"]))
ann_files = ["ann_" + os.path.basename(x) + ".txt" for x in mzML_files]


Ent = collections.namedtuple("Ent", ("Mmass name mz I adduct charge rt formu"))
Ann = collections.namedtuple(
    "Ann", ("nn dotp premz rt mz_l I_l auc feat formu m_peaks")
)

all_dat = collections.defaultdict(list)
lib_dict = dict()
for nn, ann_file in enumerate(ann_files):
    it_ann = iter(line.rstrip("\n") for line in open(ann_file))
    for line in it_ann:
        if line.startswith("NAME:"):
            name_id = [next(it_ann)]
            for line in it_ann:
                if not line.startswith("ADDUCT: "):
                    name_id.append(line)
                else:
                    name_id = (tuple(name_id), line.split(": ")[1])
                    break
            line = next(it_ann)
            premz = line[line.find(": ") + 2 :].split(", ")
            feat = True if re.fullmatch("\d+\.\d+", premz[1]) else False
            premz = (
                float(premz[1])
                if re.fullmatch("\d+\.\d+", premz[1])
                else float(premz[0])
            )
            line = next(it_ann)
            formu = line.split(": ")[1]
            line = next(it_ann)
            rt = line[line.find(": ") + 2 :].split(", ")
            rt = float(rt[1]) if re.fullmatch("\d+\.\d+", rt[1]) else float(rt[0])
            line = next(it_ann)  # peak_area
            auc = line.split(": ")[1]
            auc = float(auc) if re.fullmatch("\d+\.\d+", auc) else None
            line = next(it_ann)  # dotp
            dotp = line.split(": ")[1]
            dotp = float(dotp) if re.fullmatch("\d+\.\d+", dotp) else None
            line = next(it_ann)  # matching peaks
            m_peaks = int(line.split(": ")[1])
            next(it_ann)  # experimental_spec
            mz_list = []
            I_list = []
            for line in it_ann:
                if line and not line.startswith("LIBRARY_SPECTRUM:"):
                    mz, I = line.split()
                    mz_list.append(float(mz))
                    I_list.append(float(I))
                else:
                    break
            lib_dat = ""
            for line in it_ann:
                if line:
                    lib_dat += line + "\n"
                else:
                    break
            mz_list = [float(x) for x in mz_list]
            I_list = [x / max(I_list) * 999.0 for x in I_list]
            all_dat[name_id].append(
                Ann(nn, dotp, premz, rt, mz_list, I_list, auc, feat, formu, m_peaks)
            )
            lib_dict[name_id] = (dotp, lib_dat)

isf_dict = dict()
for name_id in all_dat:
    (name0, *_), adduct = name_id
    if name0.startswith("ISF of"):
        pmz, rt, mz = [float(x) for x in re.findall("\d+\.*\d*", name0)[:3]]
        for pmz0, rt0, mz0, adduct0 in isf_dict.keys():
            if (
                abs(pmz - pmz0) < bound_ppm(pmz * ms1ppm)
                and abs(rt - rt0) < RT_shift
                and abs(mz - mz0) < bound_ppm(mz * ms1ppm)
                and adduct == adduct0
            ):
                isf_dict[(pmz0, rt0, mz0, adduct)].append(name_id)
                break
        else:
            isf_dict[(pmz, rt, mz, adduct)] = [name_id]


name_dict0 = dict()
for key, vv in isf_dict.items():
    if len(vv) > 1:
        pmz_rt_mz = [
            tuple(float(x) for x in re.findall("\d+\.\d+", v[0])[:3]) for v, _ in vv
        ]
        key0 = statistics.median(x for x, _, _ in pmz_rt_mz)
        key1 = statistics.median(x for _, x, _ in pmz_rt_mz)
        key2 = statistics.median(x for _, _, x in pmz_rt_mz)
        isf_name = "ISF of (m/z={}, rt={}s) {}".format(key0, key1, key2)
        for v in vv:
            name_dict0[v] = ((isf_name,), key[3])


all_dat0 = collections.defaultdict(list)
lib_dict0 = collections.defaultdict(list)
for name, dat in all_dat.items():
    all_dat0[name_dict0.get(name, name)].extend(dat)
    lib_dict0[name_dict0.get(name, name)].append(lib_dict[name])


p_mz_rt = []
for ((name, *_), _), dat in all_dat0.items():
    if name.startswith("ISF of"):
        pmz, rt = [float(x) for x in re.findall("\d+\.\d+", name)[:2]]
    else:
        pmz = dat[0][2]
        rt = statistics.median(d.rt for d in dat)
    pos0 = bisect_left(p_mz_rt, (pmz - 0.01,))
    pos1 = bisect_left(p_mz_rt, (pmz + 0.01,), lo=pos0)
    for x in p_mz_rt[pos0:pos1]:
        if abs(rt - x[1]) < RT_shift:
            break
    else:
        bisect.insort(p_mz_rt, (pmz, rt))

mz_rt_name = collections.defaultdict(list)
for name, dat in all_dat0.items():
    if name[0][0].startswith("ISF of"):
        pmz, rt = [float(x) for x in re.findall("\d+\.\d+", name[0][0])[:2]]
    else:
        pmz = dat[0][2]
        rt = statistics.median(d.rt for d in dat)
    pos0 = bisect_left(p_mz_rt, (pmz - 0.01,))
    pos1 = bisect_left(p_mz_rt, (pmz + 0.01,), lo=pos0)
    for n, x in enumerate(p_mz_rt[pos0:pos1], pos0):
        if abs(rt - x[1]) < RT_shift:
            mz_rt_name[n].append(name)
            break
    else:
        print("err")
        sys.exit()

for key, names in sorted(mz_rt_name.items()):
    if all(x[0][0].startswith("ISF of ") for x in names):
        del mz_rt_name[key]


def eligible_parent(x):
    return not x[0][0].startswith("ISF of ")


for key, names in sorted(mz_rt_name.items()):
    if any(x[0][0].startswith("ISF of ") for x in names) and all(
        not eligible_parent(x) for x in names
    ):
        mz_rt_name[key] = [x for x in names if not x[0][0].startswith("ISF of ")]

isf_mz_rt = []
for key, names in sorted(mz_rt_name.items()):
    for (name, *_), _ in names:
        if name.startswith("ISF of"):
            pmz, rt, mz = [float(x) for x in re.findall("\d+\.\d+", name)[:3]]
            isf_mz_rt.append((mz, rt, pmz))
isf_mz_rt.sort()


grp_count = 1
name_mz_rt = dict()
name_isf_grp = dict()
possible_isf = set()
for key, names in sorted(mz_rt_name.items()):
    flag = False
    for name in (x for x in names if not x[0][0].startswith("ISF of ")):
        if name not in all_dat0:
            print("daf")
        mz0 = statistics.median(x[2] for x in all_dat0[name])
        rt0 = statistics.median(x[3] for x in all_dat0[name])
        pos0 = bisect_left(isf_mz_rt, (mz0 - bound_ppm(mz0 * ms1ppm),))
        pos1 = bisect_left(isf_mz_rt, (mz0 + bound_ppm(mz0 * ms1ppm),))
        for mz, rt, _ in isf_mz_rt[pos0:pos1]:
            if abs(rt - rt0) < RT_shift:
                flag = True
                possible_isf.add(name)
                break
    isf_in_names = [x for x in names if x[0][0].startswith("ISF of")]
    if flag and len(isf_in_names) > -1:
        continue

    for name in names:
        name_mz_rt[name] = grp_count
    grp_count += 1
    isf_pres = any(x[0][0].startswith("ISF of") for x in names)
    for name in names:
        name_isf_grp[name] = isf_pres


mz_rt_n_dict = {
    k: collections.defaultdict(list) for k in list(adduct_list)
}  # +['isf']}
mz_rt_n = []
for k, n in name_mz_rt.items():
    if k[0][0].startswith("ISF of"):
        mz_rt_n.append((n, [k]))
    else:
        mz_rt_n_dict[list(mz_rt_n_dict)[0]][n].append(k)


n0_names = dict()
for n0, names in sorted(
    [
        inner
        for outer in [list(x.items()) for x in mz_rt_n_dict.values()]
        for inner in outer
    ]
    + mz_rt_n
):
    if not names[0][0][0].startswith("ISF of "):
        n0_names[n0] = [x[0][0] for x in names]


isf_dat = set()
for n0, names in sorted(
    [
        inner
        for outer in [list(x.items()) for x in mz_rt_n_dict.values()]
        for inner in outer
    ]
    + mz_rt_n
):
    if names[0][0][0].startswith("ISF of "):
        for d in all_dat0[names[0]]:
            isf_dat.add((d.nn, d.premz, d.rt))


with open("ann_All.txt", "w") as cpd_ann, open("quant_All.txt", "w") as quant_auc, open(
    "All.txt", "w"
) as cpd_una:
    quant_auc.write(
        "group\tISF\tname\tadduct\tformula\tfeature_m/z\tmatching_peaks(median)\tMin.\t1st Qu.\tMedian\t3rd Qu.\tMax.\t%detected\t"
    )
    quant_auc.write("\t".join(x[:-5] for x in mzML_files) + "\t")
    quant_auc.write("\t".join("RT_" + x[:-5] for x in mzML_files) + "\t")
    quant_auc.write("\t".join("score_" + x[:-5] for x in mzML_files) + "\n")
    n1 = 1
    prev_n0 = 1
    firstflag = False
    for n0, names in sorted(
        [
            inner
            for outer in [list(x.items()) for x in mz_rt_n_dict.values()]
            for inner in outer
        ]
        + mz_rt_n
    ):
        dat = []
        dat0 = collections.defaultdict(list)
        lib_dict1 = []
        for name in names:
            auc_dfdict = collections.defaultdict(list)
            for d in all_dat0[name]:
                if (
                    name[0][0].startswith("ISF of ")
                    or (d.nn, d.premz, d.rt) not in isf_dat
                ):
                    auc_dfdict[d.nn].append(d)

            lib_dict1.extend([x + (name,) for x in lib_dict0[name]])
            for nn in range(len(mzML_files)):
                dotp_rt_auc = auc_dfdict.get(nn)

                if dotp_rt_auc:
                    for x in dotp_rt_auc:
                        dat0[nn].append(x)

        for nn in dat0.keys():
            if any(x.feat for x in dat0[nn]):
                dat0[nn] = [x for x in dat0[nn] if x.feat]
            else:
                dat0[nn] = [max(dat0[nn], key=operator.attrgetter("dotp"))]

        if not dat0 or (
            sum(x.feat for xx in dat0.values() for x in xx) == 0
            and len(dat0) < len(mzML_files) * 0.8
        ):
            continue
        for nn, dotp_l in dat0.items():
            sdotp = sorted(dotp_l, key=operator.attrgetter("dotp"), reverse=True)
            dat.append(sdotp[0])
            for p in sdotp[1:]:
                if sdotp[0].dotp - p.dotp < 0.1:
                    dat.append(p)
                else:
                    break
        dat.sort()

        nameset = sorted(set(x for name in names for x in name[0]))
        if any((name in possible_isf) for name in names):
            nameset = ["possibly an ISF"] + nameset
        adductset = sorted(set(name[1] for name in names))

        premz_l = [d.premz for d in dat]
        rt_l = [d.rt for d in dat]
        cpd_ann.write("NAME:\n")
        cpd_ann.write(max(lib_dict1)[2][0][0] + "\n")
        cpd_ann.write("ADDUCT: " + max(lib_dict1)[2][1] + "\n")
        cpd_ann.write("SAMPLE, RT, DOT_PRODUCT, PEAK_AREA\n")
        for d in dat:
            cpd_ann.write(
                "{:52}{:<9.2f}{:<6.2f}{}\n".format(
                    mzML_files[d.nn][-55:-5],
                    d.rt,
                    d.dotp,
                    format(d.auc, ".1f") if d.auc else "no_ms1_feature",
                )
            )
        cpd_ann.write("PRECURSOR_M/Z: {:.5f}\n".format(statistics.median(premz_l)))
        cpd_ann.write(
            "FORMULA: {}\n".format(",".join(sorted(set(d.formu for d in dat))))
        )
        cpd_ann.write("RT: {:.2f}\n".format(statistics.median(rt_l)))
        cpd_ann.write("EXPERIMENTAL_SPECTRUM:\n")
        max_dat = max(dat, key=operator.attrgetter("dotp"))
        for mz, I in zip(max_dat[4], max_dat[5]):
            cpd_ann.write("{:.5f} {:.2f}\n".format(mz, I))
        cpd_ann.write("LIBRARY_SPECTRUM:\n")
        lib_dat = max(lib_dict1)[1]
        cpd_ann.write(lib_dat)
        cpd_ann.write("\n")

        cpd_una.write("NAME: {}\n".format(max(lib_dict1)[2][0][0]))
        cpd_una.write("PRECURSORTYPE: [{}]+\n".format(max(lib_dict1)[2][1]))
        cpd_una.write("PRECURSORMZ: {:.5f}\n".format(statistics.median(premz_l)))
        cpd_una.write("RETENTIONTIME: {:.2f}\n".format(statistics.median(rt_l) / 60))
        cpd_una.write("Num Peaks: {}\n".format(min(len(max_dat[5]), 20)))
        for mz, I in zip(max_dat[4], max_dat[5][:20]):
            cpd_una.write("{:.5f} {:.2f}\n".format(mz, I))
        cpd_una.write("\n")

        if prev_n0 != n0 and firstflag:
            n1 += 1
        firstflag = True
        quant_auc.write(str(n1) + "\t")
        quant_auc.write(("*" if name_isf_grp[name] else "") + "\t")
        if nameset[0].startswith("ISF of "):
            quant_auc.write("ISF of " + " --- ".join(n0_names[n0]))
        else:
            quant_auc.write(" --- ".join(nameset))
        prev_n0 = n0
        quant_auc.write(
            "\t{}\t{}\t{:.5f}\t{:.1f}".format(
                ",".join(adductset),
                ",".join(sorted(set(d.formu for d in dat))),
                statistics.median(premz_l),
                statistics.median(d.m_peaks for d in dat),
            )
        )

        s_rt = sorted(dat_n.rt for dat_n in dat)
        Qu1 = s_rt[round((len(s_rt) - 1) * 0.25)] / 60
        Qu2 = s_rt[round((len(s_rt) - 1) * 0.50)] / 60
        Qu3 = s_rt[round((len(s_rt) - 1) * 0.75)] / 60
        quant_auc.write(
            "\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                s_rt[0] / 60,
                Qu1,
                Qu2,
                Qu3,
                s_rt[-1] / 60,
                len({x.nn for x in dat}) / len(mzML_files),
            )
        )

        line_str = []
        for nn in range(len(mzML_files)):
            pos0 = bisect_left([x.nn for x in dat], nn)
            pos1 = bisect.bisect_right([x.nn for x in dat], nn, lo=pos0)
            line_str.append(
                (
                    ",".join(
                        format(dat_n.auc, ".1f") if dat_n.auc else "no_ms1_feature"
                        for dat_n in dat[pos0:pos1]
                    ),
                    ",".join(format(dat_n.rt / 60, ".2f") for dat_n in dat[pos0:pos1]),
                    ",".join(format(dat_n.dotp, ".2f") for dat_n in dat[pos0:pos1]),
                )
            )
        quant_auc.write(
            "\t"
            + "\t".join(x for x, _, _ in line_str)
            + "\t"
            + "\t".join(x for _, x, _ in line_str)
            + "\t"
            + "\t".join(x for _, _, x in line_str)
            + "\n"
        )

    una_files = ["una_" + os.path.basename(x) + ".txt" for x in mzML_files]
    una_dict = dict()
    for ii, una_file in enumerate(una_files):
        it_cpd = iter(line.rstrip() for line in open(una_file))
        for line in it_cpd:
            lsp = line.split(": ", 1)
            if lsp[0] == "PRECURSORMZ":
                ms1mz = float(lsp[1])
            elif lsp[0] == "RETENTIONTIME":
                RT = float(lsp[1]) * 60
            elif lsp[0] == "Num Peaks":
                frag_mz = []
                frag_I = []
                for lsp in (line.split() for line in it_cpd):
                    if not lsp:
                        break
                    frag_mz.append(float(lsp[0]))
                    frag_I.append(float(lsp[1]))
                una_dict[(ms1mz, RT, ii)] = (frag_mz, frag_I)

    s_mz_rt = sorted(una_dict)
    i_mz_rt = set()
    grp_mz_rt = []
    for pos0, (pmz, rt, ii) in enumerate(s_mz_rt):
        if (pmz, rt, ii) not in i_mz_rt:
            mz_rt_l = []
            pos1 = bisect_left(s_mz_rt, (pmz + 0.01,), lo=pos0)
            for x in s_mz_rt[pos0:pos1]:
                if x not in i_mz_rt and abs(rt - x[1]) < RT_shift:
                    i_mz_rt.add(x)
                    mz_rt_l.append(x)
            grp_mz_rt.append(mz_rt_l)
    grp_mz_rt = [grp for grp in grp_mz_rt if len(grp) > len(una_files) / 2]
    lib_ent = DDAreadlib.get_cpds()
    for grp in grp_mz_rt:
        max_spec = max([(x, una_dict[x]) for x in grp], key=lambda x_y: max(x_y[1][1]))
        median_mz = statistics.median(x[0] for x in grp)
        median_rt = statistics.median(x[1] for x in grp) / 60

        err_bd = bound_ppm(median_mz * ms1ppm)
        pos_0 = bisect_left(lib_ent, (median_mz - err_bd,))
        pos_1 = bisect_left(lib_ent, (median_mz + err_bd,), lo=pos_0)
        lib_ent_ = lib_ent[pos_0:pos_1]
        if lib_ent_:
            s_lib_ent_ = sorted(lib_ent_, key=lambda x: abs(x.Mmass - median_mz))
            cpd_una.write(
                "NAME: (putative) {}\n".format(
                    " --- ".join(list({x.name: None for x in s_lib_ent_}))
                )
            )
        else:
            cpd_una.write("NAME: {:.4f}_{:.2f}\n".format(median_mz, median_rt))
        cpd_una.write("PRECURSORTYPE: [unknown]+\n")
        cpd_una.write("PRECURSORMZ: {:.5f}\n".format(median_mz))
        cpd_una.write("RETENTIONTIME: {:.2f}\n".format(median_rt))
        cpd_una.write("Num Peaks: {}\n".format(len(max_spec[1][1])))
        for mz, i in zip(*max_spec[1]):
            cpd_una.write("{:.5f} {:.2f}\n".format(mz, i))
        cpd_una.write("\n")

print("done")
