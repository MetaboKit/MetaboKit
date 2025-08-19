import collections
import os
import sys
import re
import glob

import DDAcommonfn

Ent = collections.namedtuple("Ent", ("Mmass name mz I adduct charge rt formu"))
param_set = {
    "mzML_files",
    "library",
}
param_dict = DDAcommonfn.read_param(param_set)

lib_types = [x.split(" #", 1)[0].strip() for x in param_dict["library"].splitlines()]

mzML_files = sorted(glob.glob(param_dict["mzML_files"]))

for line in open(mzML_files[0]):
    if "MS:1000130" in line:
        ispos = True
        break
    elif "MS:1000129" in line:
        ispos = False
        break

libpaths = []
script_dir = os.path.abspath(os.path.dirname(sys.argv[0])) + "/libs/"
if "nist" in lib_types:
    libpaths.append(script_dir + "nist.Original.bak.msp")
if "LipidBlast" in lib_types:
    if ispos:
        libpaths.append(script_dir + "LipidBlast-ASCII-spectra/LipidBlast-pos.msp")
        libpaths.extend(
            glob.glob(script_dir + "LipidBlast-ASCII-spectra/custom-libs/*pos.msp")
        )
    else:
        libpaths.append(script_dir + "LipidBlast-ASCII-spectra/LipidBlast-neg.msp")
        libpaths.extend(
            glob.glob(script_dir + "LipidBlast-ASCII-spectra/custom-libs/*neg.msp")
        )
if "LipidBlast-fork" in lib_types:
    if ispos:
        libpaths.append(script_dir + "MSDIAL-TandemMassSpectralAtlas-VS68-Pos.msp")
    else:
        libpaths.append(script_dir + "MSDIAL-TandemMassSpectralAtlas-VS68-Neg.msp")
if "hmdb" in lib_types:
    libpaths.append(script_dir + "hmdb_metabolites.xml")
if "msdial" in lib_types:
    if ispos:
        libpaths.append(script_dir + "MSMS-Public-Pos-VS15.msp")
    else:
        libpaths.append(script_dir + "MSMS-Public-Neg-VS15.msp")
if "sling" in lib_types:
    if ispos:
        libpaths.append(script_dir + "Accurate_Mass_MRM_list_v1.txt")
if "NoMatch" in lib_types[0]:
    libpaths = [script_dir + "Database_Dec2017.txt"]
    if len(lib_types[0]) > 9:
        libpaths = [lib_types[0].split(" ", 1)[1]]
        print(lib_types, libpaths)

metabokit = []
if any(x.startswith("user ") for x in lib_types):
    for x in lib_types:
        if x.startswith("user "):
            libpaths.append(x)


cpd_list = {
    "CO": 27.99491,
    "H": 1.007825,
    "Li": 7.016004,
    "NH4": 18.03437,
    "Na": 22.98977,
    "Na2": 45.97954,
    "H": 1.007825,
    "2H": 2.015650,
    "2H": 2.015650,
    "2I": 253.8089,
    "2K": 77.92741,
    "2Na": 45.97954,
    "3H": 3.023475,
    "3K": 116.8911,
    "3Na": 68.96931,
    "H": 1.007825,
    "H2O": 18.01056,
    "I": 126.9045,
    "K": 38.96371,
    "NH3": 17.02655,
    "NH4": 18.03437,
    "Na": 22.98977,
    "OH": 17.00274,
    "": 0,
    "FA": 46.00548,
    "Hac": 60.02113,
    "ACN": 41.02655,
}


def read_lib(libpath):
    adductset = set()
    lib_dict = collections.defaultdict(list)
    libpath0 = libpath.split("/libs/")[-1]
    if "LipidCreatorValidStudy_MRM_Workklist_V1_forHyungWon.txt" in libpath0:
        qq = open(libpath)
        qq.readline()
        nnnn = 0
        for nnn, line in enumerate(qq):
            lsp = line.split("\t")
            name = "SLING_MRM_list " + lsp[1] + " " + lsp[11]
            ms1mz = lsp[2]
            adduct = lsp[11]
            charge = adduct[adduct.rfind("]") + 1]
            charge = charge if charge.isdigit() else "1"
            adduct = adduct[adduct.find("[") + 1 : adduct.rfind("]")].replace(" ", "")
            adductset.add(adduct)
            adduct = tuple(x.strip() for x in re.split("(\-|\+)", adduct))

            frag_mz = [lsp[3]]
            frag_I = ["1"]
            if frag_mz:
                lib_dict[
                    (ms1mz, charge, tuple(frag_mz), tuple(frag_I), adduct, "NA", "")
                ].append(name)
                nnnn += 1
        print(list(lib_dict)[0])
        print(adductset)
        print(nnn, nnnn)
    if "Accurate_Mass_MRM_list_v1.txt" in libpath0:
        slingi = open(libpath)
        next(slingi)
        for line in slingi:
            lsp = line.rstrip("\n").split("\t")
            ms1mz = lsp[5]
            iont = lsp[4]
            adduct = iont[iont.find("[M") + 1 : iont.rfind("]")]
            adduct = tuple(re.split("(\-|\+)", adduct))
            charge = 1
            frag_mz = (lsp[6],)
            frag_I = ("1",)
            formu = lsp[3]
            name = lsp[2] + " SLING_MRM"
            lib_dict[(ms1mz, charge, frag_mz, frag_I, adduct, "NA", formu)].append(name)
    if "LipidBlast-ASCII-spectra" in libpath0:
        it_cpd = iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            lsp = line.split(": ", 1)
            if lsp[0] == "Name":
                name = lsp[1] + " LipidBlast"
                adduct = line[line.find("[M") + 1 : line.rfind("]")]
                adduct = tuple(re.split("(\-|\+)", adduct))
                charge = line[line.rfind("]") + 1 :]
                charge = charge[1] if charge[0] == "(" else charge[0]
                charge = charge if charge.isdigit() else "1"
            elif lsp[0] == "PRECURSORMZ":
                ms1mz = lsp[1]
            elif lsp[0] == "Comment":
                formu = lsp[1].rsplit("; ", 1)[-1]
            elif lsp[0] == "Num Peaks":
                frag_mz = []
                frag_I = []
                for lsp in (line.split() for line in it_cpd):
                    if not lsp:
                        break
                    frag_mz.append(lsp[0])
                    frag_I.append(lsp[1])
                lib_dict[
                    (ms1mz, charge, tuple(frag_mz), tuple(frag_I), adduct, "NA", formu)
                ].append(name)

    if "nist.Original.bak.msp" in libpath0:
        it_cpd = iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            lsp = line.split(": ", 1)
            if lsp[0] == "Name":
                name = lsp[1]
                line = next(it_cpd)
                name += " " + line[line.find("[") :] + " NIST"
                adduct = line[line.find("[") + 1 : line.rfind("]")]
                adduct = tuple(re.split("(\-|\+)", adduct))
                charge = line[line.rfind("]") + 1]
                charge = charge if charge.isdigit() else "1"
                posneg = line[-1]
            elif lsp[0] == "Formula":
                formu = lsp[1]
            elif lsp[0] == "PrecursorMZ":
                ms1mz = lsp[1]
            elif lsp[0] == "Num peaks":
                frag_mz = []
                frag_I = []
                for lsp in (line.split() for line in it_cpd):
                    if not lsp:
                        break
                    frag_mz.append(lsp[0])
                    frag_I.append(lsp[1])
                if re.fullmatch("\d+\.\d{2,}", frag_mz[0]) is None:
                    continue  # remove if <n dp
                if (ispos and posneg != "-") or (not ispos and posneg == "-"):
                    lib_dict[
                        (
                            ms1mz,
                            charge,
                            tuple(frag_mz),
                            tuple(frag_I),
                            adduct,
                            "NA",
                            formu,
                        )
                    ].append(name)
    if (
        "MSDIAL-TandemMassSpectralAtlas-" in libpath0
        or "MSMS-Public-" in libpath0
        or libpath0.startswith("user ")
    ):
        lib = ""
        formu = ""
        if "MSDIAL-TandemMassSpectralAtlas-" in libpath0:
            lib = " Atlas"
        elif "MSMS-Public-" in libpath0:
            lib = " MSDIAL"
        if libpath0.startswith("user "):
            libpath = libpath0.split(" ", 1)[1]
            if not open(libpath):
                print(libpath + " not found")
                sys.exit()
        it_cpd = iter(line.rstrip() for line in open(libpath))
        for line in it_cpd:
            lsp = line.split(": ", 1)
            if lsp[0] == "NAME":
                name = lsp[1]
                RT = "NA"
            elif lsp[0] == "PRECURSORMZ":
                ms1mz = lsp[1]
            elif line.startswith("PRECURSORTYPE: "):
                adduct = line[line.find("[") + 1 : line.rfind("]")]
                name += line.split(":", 1)[1] + lib
                adduct = tuple(re.split("(\-|\+)", adduct))
                pos = line.rfind("]")
                if pos != -1:
                    charge = line[pos + 1]
                    charge = "1" if any(charge == ch for ch in "+-") else charge
                else:
                    charge = "1"
            elif lsp[0] == "FORMULA":
                formu = lsp[1]
            elif lsp[0] == "RETENTIONTIME" and libpath0.startswith("user "):
                RT = str(float(lsp[1]) * 60)
            elif lsp[0] == "Num Peaks":
                frag_mz = []
                frag_I = []
                for lsp in (line.split() for line in it_cpd):
                    if not lsp:
                        break
                    frag_mz.append(lsp[0])
                    frag_I.append(lsp[1])
                if frag_mz:
                    lib_dict[
                        (
                            ms1mz,
                            charge,
                            tuple(frag_mz),
                            tuple(frag_I),
                            adduct,
                            RT,
                            formu,
                        )
                    ].append(name)
    print(libpath0)
    print(len(lib_dict))
    return lib_dict


def get_cpds():  # from libs
    lib_dict = dict()
    for libpath in libpaths:
        lib_dict.update(read_lib(libpath))
    lib_ent = []
    for ent, name in lib_dict.items():
        ms1mz, charge, frag_mz, frag_I, adduct, RT, formu = ent
        ms1mz, charge = float(ms1mz), int(charge)
        RT = None if RT == "NA" else float(RT)
        frag_mz = [float(x) for x in frag_mz]
        frag_I = [float(x) for x in frag_I]

        I_mz_list = sorted(zip(frag_I, frag_mz), reverse=True)[:10]
        frag_mz = tuple(x for _, x in I_mz_list)
        frag_I = tuple(x for x, _ in I_mz_list)

        name = ", ".join(sorted(set(name)))

        lib_ent.append(
            Ent(ms1mz, name, frag_mz, frag_I, "".join(adduct), charge, RT, formu)
        )
        continue

        Mmass = ms1mz * charge
        adduct = adduct.split(",")
        for pm, cpd in zip(adduct[1::2], adduct[2::2]):
            if pm == "+":
                Mmass -= cpd_list[cpd]
            elif pm == "-":
                Mmass += cpd_list[cpd]
            else:
                sys.exit()
        if adduct[0][0].isdigit():  # M,2M,3M
            Mmass /= int(adduct[0][0])

        lib_ent.append(Ent(Mmass, name, frag_mz, frag_I, None, None, RT, formu))
    return sorted(lib_ent)
