
import pickle
import requests
import json
import re
import os

#from https://www.genome.jp/kegg-bin/show_brite?ko00002.keg
kegg_module='modules.json'
pickle_path='modules.pickle'

def save_metrics(pickle_path,to_pickle):
    with open(pickle_path, 'wb') as handle:
        pickle.dump(to_pickle, handle)

def load_metrics(pickle_path):
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as handle:
            pickled_results= pickle.load(handle)
            return pickled_results

def find_ko(string_to_search):
    res=set()
    #I could do upper and lower case but since it's only one letter, it's not very safe...
    pattern = re.compile('(?<![A-Za-z])K\d{4,}')
    search= re.finditer(pattern,string_to_search)
    for i in search:
        res.add(i.group())
    return res

def get_ko_from_module(module_id):
    url = 'https://www.genome.jp/dbget-bin/www_bget?md:' + module_id
    webpage = None
    c = 0
    while not webpage and c <= 10:
        req = requests.get(url)
        try:
            webpage = req.text
        except:
            c += 1
        start=re.search('<nobr>Definition</nobr>',webpage).span()[1]
        webpage=webpage[start:]
        end=re.search('</div></div></td></tr>',webpage).span()[0]
        webpage=webpage[:end]
        ko_set=find_ko(webpage)
    return ko_set



def read_modules(file_path):
    with open(file_path) as file: json_modules= json.load(file)['children'][0]['children']
    tree_modules={}
    for main_path in json_modules:
        main_path_name=main_path['name']
        if main_path_name not in tree_modules: tree_modules[main_path_name]={}
        sub_pathways=main_path['children']
        for sub_path in sub_pathways:
            sub_path_name=sub_path['name']
            if sub_path_name not in tree_modules[main_path_name]: tree_modules[main_path_name][sub_path_name] = {}
            modules=sub_path['children']
            for module in modules:
                module_name=module['name'].split()[0]
                tree_modules[main_path_name][sub_path_name][module_name]=get_ko_from_module(module_name)
    return tree_modules