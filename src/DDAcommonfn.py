import collections


def read_param(param_set):
    param_dict = dict()
    with open("param.txt") as input_param:
        key_ = ""
        value_ = ""
        for line in (l.strip() for l in input_param if l[0] != "#"):
            if not line and key_ and value_:
                param_dict[key_] = value_[:-1]
                key_ = ""
                value_ = ""
            elif line in param_set:
                key_ = line
            elif key_:
                value_ += line + "\n"
        param_dict[key_] = value_[:-1]
    return param_dict


def bound_ppm(mz_diff):
    return max(0.003, mz_diff)


Ent = collections.namedtuple("Ent", ("Mmass name adduct rt"))


def read_msp(libpath):
    lib_l = []
    lib_dict = collections.defaultdict(list)
    libpath0 = libpath.split("/libs/")[-1]
    if libpath0.startswith("metabokit "):
        libpath = libpath0.split(" ", 1)[1]
        if not open(libpath):
            print(libpath + " not found")
            sys.exit()
    it_cpd = iter(line.rstrip() for line in open(libpath))
    for line in it_cpd:
        if line.upper().startswith("NAME: "):
            name = line.split(": ", 1)[1]
            RT = "NA"
            if "[M" in line:
                adduct = line[line.find("[M") + 1 : line.rfind("]")]
            else:
                adduct = ""
        elif line.upper().startswith("PRECURSORMZ: "):
            ms1mz = float(line.split(": ", 1)[1])
        elif line.startswith("PRECURSORTYPE: "):
            adduct = line.split(": ", 1)[1]
            if "[M" in line:
                adduct = line[line.find("[M") + 1 : line.rfind("]")]
            else:
                adduct = ""
        elif line.startswith("RETENTIONTIME: ") and libpath0.startswith("metabokit "):
            RT = float(line.split(": ", 1)[1]) * 60
        elif not line:
            lib_dict[(ms1mz, adduct, RT)].append(name)
    if (ms1mz, adduct, RT) not in lib_dict or name not in lib_dict[(ms1mz, adduct, RT)]:
        lib_dict[(ms1mz, adduct, RT)].append(name)
    for (ms1mz, adduct, RT), name in lib_dict.items():
        lib_l.append(Ent(ms1mz, "---".join(name), adduct, RT))
    print(libpath, len(lib_l))
    return lib_l
