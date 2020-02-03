


def read_param(param_set):
    param_dict=dict()
    with open('param.txt') as input_param:
        key_=''
        value_=''
        for line in (l.strip() for l in input_param if l[0]!='#'):
            if not line and key_ and value_:
                param_dict[key_]=value_[:-1]
                key_=''
                value_=''
            elif line in param_set:
                key_=line
            elif key_:
                value_+=line+'\n'
    return param_dict


def bound_ppm(mz_diff):
    return max(.003,mz_diff)

