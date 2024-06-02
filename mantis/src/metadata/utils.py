import re


def find_tcdb(string_to_search):
    res = set()
    tc_pattern = re.compile(r'(?<![A-Za-z])\(TC\s\d\.[A-Z\-](\.(\d+|\-)){1,2}\)')
    search = re.finditer(tc_pattern, string_to_search)
    for i in search:
        res.add(i.group())
    return res


def find_ko(string_to_search):
    res = set()
    # I could do upper and lower case but since it's only one letter, it's not very safe...
    pattern = re.compile(r'(?<![A-Za-z])K\d{4,}')
    search = re.finditer(pattern, string_to_search)
    for i in search:
        res.add(i.group())
    return res


def find_pfam(string_to_search):
    res = set()
    duf = re.compile(r'(?<![A-Za-z])(DUF|duf)\d{2,}')
    pfam = re.compile(r'(?<![A-Za-z])((U|u)?PF|pf)\d{3,}')
    for pattern in [duf, pfam]:
        search = re.finditer(pattern, string_to_search)
        for i in search:
            res.add(i.group())
    return res


def find_cog(string_to_search):
    res = set()
    pattern = re.compile(r'(?<![A-Za-z])(COG|cog)\d{3,}')
    search = re.finditer(pattern, string_to_search)
    for i in search:
        res.add(i.group())
    return res


def find_arcog(string_to_search):
    res = set()
    pattern = re.compile(r'(AR|ar)(COG|cog)\d{3,}')
    search = re.finditer(pattern, string_to_search)
    for i in search:
        res.add(i.group())
    return res


def find_tigrfam(string_to_search):
    res = set()
    pattern = re.compile(r'(?<![A-Za-z])(TIGR)\d{3,}')
    search = re.finditer(pattern, string_to_search)
    for i in search:
        res.add(i.group())
    return res


def find_go(string_to_search):
    res = set()
    pattern = re.compile(r'(?<![A-Za-z])GO\d{3,}')
    search = re.finditer(pattern, string_to_search)
    for i in search:
        res.add(i.group())
    return res



def find_ecs(string_to_search, required_level=3):
    res = set()
    # greedy match of confounders
    ec_pattern = re.compile(r'\d(\.(-|\d{1,3}|([a-zA-Z]\d{1,3}))){2,3}')
    search = re.finditer(ec_pattern, string_to_search)
    for i in search:
        ec = i.group()
        passed = False
        start = i.span()[0]
        end = i.span()[1]
        if len(string_to_search) > end:
            if string_to_search[start - 1] != '.' and \
                    string_to_search[end] != '.' \
                    and not re.match(r'\.|[a-zA-Z]|\d{1,3}', string_to_search[end]) and not re.match(r'-',string_to_search[end]):
                passed = True
        else:
            if string_to_search[start - 1] != '.':
                passed = True
        if passed:
            if ec.count('.') >= required_level - 1:
                if ec.count('.') + 1 - ec.count('-') >= required_level:
                    res.add(ec)
    return res

def get_common_links_metadata(input_string: str='', metadata_dict: dict={}) -> dict:
    if not input_string:
        return metadata_dict
    metadata_function_mapping = {
        'enzyme_ec': find_ecs,
        'tcdb': find_tcdb,
        'tigrfam': find_tigrfam,
        'kegg_ko': find_ko,
        'pfam': find_pfam,
        'cog': find_cog,
        'arcog': find_arcog,
        'go': find_go,
        }
    for metadata_str, metadata_function in metadata_function_mapping.items():
        metadata_set = metadata_function(input_string)
        if metadata_str not in metadata_dict:
            metadata_dict[metadata_str] = set()
        metadata_dict[metadata_str].update(metadata_set)
