'''
This file contains helper functions neededby the manipulator
They are in a separate file such that they can also be used separately
=> esp. bcs then Kratos does not have to be imported
'''


def DictToPrettyString(o, level=0):
    '''
    This function returns a nicely formatted string of a dictionary
    This can e.g. be written to a file
    '''
    # source: https://stackoverflow.com/questions/10097477/python-json-array-newlines
    INDENT = 4
    SPACE = " "
    NEWLINE = "\n"
    ret = ""
    if isinstance(o, dict):
        if len(o) == 0:
            ret += "{}"
        elif len(o) == 1: # and list(o.keys())[0] == "@table"
            ret += "{"
            for k in o.keys():
                v = o[k]
                ret += '"' + str(k) + '" : '
                ret += DictToPrettyString(v, level + 1)

            ret += "}"
        else:
            ret += "{" + NEWLINE
            comma = ""
            for k in reversed(sorted(o.keys())): # reversed to make i better readable
                v = o[k]
                ret += comma
                comma = ",\n"
                ret += SPACE * INDENT * (level+1)
                ret += '"' + str(k) + '" :' + SPACE
                ret += DictToPrettyString(v, level + 1)

            ret += NEWLINE + SPACE * INDENT * level + "}"
    elif isinstance(o, str):
        ret += '"' + o + '"'
    elif isinstance(o, list):
        if len(o) > 0:
            if isinstance(o[0], (int,  float, str, dict)): # dict is needed for the composite-materials
                ret += NEWLINE + SPACE * INDENT * (level+1) + "["
                num_vals = len(o)
                for i in range(num_vals):
                    entry = o[i]
                    ret += DictToPrettyString(entry)
                    if i <= num_vals-2:
                        ret += ", "
                    else:
                        ret += "]"
            else:
                ret += "[" + ",".join([DictToPrettyString(e, level+1) for e in o]) + "]"
        else:
            ret += "[]"
    elif isinstance(o, bool):
        ret += "true" if o else "false"
    elif isinstance(o, int):
        ret += str(o)
    elif isinstance(o, float):
        ret += '%.7g' % o
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.integer):
        ret += "[" + ','.join(map(str, o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.inexact):
        ret += "[" + ','.join(map(lambda x: '%.7g' % x, o.flatten().tolist())) + "]"
    elif o is None:
        ret += 'null'
    elif isinstance(o, TabularData):
        err
        data = o.GetTabularData()
        ret += "[" + NEWLINE
        num_rows = len(data)
        for i in range(num_rows-1):
            row = data[i]
            ret += SPACE * INDENT * (level+1) + "[" + str(row[0]) + " , " + str(row[1]) + "]," + NEWLINE
        ret += SPACE * INDENT * (level+1) + "[" + str(row[0]) + " , " + str(row[1]) + "]" + NEWLINE
        ret += SPACE * INDENT * (level+1) + "]"
    else:
        raise TypeError("Unknown type '%s' for json serialization" % str(type(o)))
    return ret