
import copy
import json
import re
import requests
import pickle


#this class is based on https://github.com/merenlab/anvio
class KeggTreeGenerator():

    def split_path(self, step):
        """This function handles compound steps that should be split into multiple alternative paths.
        It first splits the input step into substeps, and then since each substep could be its own mini-definition,
        it recursively calls the definition unrolling function to parse it. The list of all alternative paths
        that can be made from this step is returned.
        """

        if step[0] == "(" and step[-1] == ")":
            substeps = self.split_by_delim_not_within_parens(step[1:-1], ",")
            if not substeps: # if it doesn't work, try without removing surrounding parentheses
                substeps = self.split_by_delim_not_within_parens(step, ",")
        else:
            substeps = self.split_by_delim_not_within_parens(step, ",")

        alt_path_list = []
        for s in substeps:
            alt_paths_from_substep = self.recursive_definition_unroller(s)
            for a in alt_paths_from_substep:
                alt_path_list.append(a)

        return alt_path_list

    def split_by_delim_not_within_parens(self, d, delims, return_delims=False):
        """Takes a string, and splits it on the given delimiter(s) as long as the delimeter is not within parentheses.
        This function exists because regular expressions don't handle nested parentheses very well. It is used in the
        recursive module definition unrolling functions to split module steps, but it is generically written in case
        it could have other uses in the future.
        The function can also be used to determine if the parentheses in the string are unbalanced (it will return False
        instead of the list of splits in this situation)
        PARAMETERS
        ==========
        d : str
            string to split
        delims : str or list of str
            a single delimiter, or a list of delimiters, to split on
        return_delims : boolean
            if this is true then the list of delimiters found between each split is also returned
        RETURNS
        =======
        If parentheses are unbalanced in the string, this function returns False. Otherwise:
        splits : list
            strings that were split from d
        delim_list : list
            delimiters that were found between each split (only returned if return_delims is True)
        """

        parens_level = 0
        last_split_index = 0
        splits = []
        delim_list = []
        for i in range(len(d)):
            # only split if not within parentheses
            if d[i] in delims and parens_level == 0:
                splits.append(d[last_split_index:i])
                delim_list.append(d[i])
                last_split_index = i + 1 # we add 1 here to skip the space
            elif d[i] == "(":
                parens_level += 1
            elif d[i] == ")":
                parens_level -= 1

            # if parentheses become unbalanced, return False to indicate this
            if parens_level < 0:
                return False
        splits.append(d[last_split_index:len(d)])

        if return_delims:
            return splits, delim_list
        return splits

    def recursive_definition_unroller(self,step):
        """This function recursively splits a module definition into its components.
        First, the definition is split into its component steps (separated by spaces).
        Each step is either an atomic step (a single KO, module number, '--', or nonessential KO starting with '-'),
        a protein complex, or a compound step.
        Atomic steps are used to extend each path that has been found so far. Protein complexes are split into
        their respective components, which may be split further by the split_paths() function to find all possible
        alternative complexes, before being used to extend each path. Compound steps are split and recursively processed
        by the split_paths() function before the resulting downstream paths are used to extend each path.
        PARAMETERS
        ==========
        step : str
            step definition to split into component steps as necessary
        RETURNS
        =======
        paths_list : list
            all paths that the input step has been unrolled into
        """

        split_steps = self.split_by_delim_not_within_parens(step, " ")
        paths_list = [[]]  # list to save all paths, with initial empty path list to extend from
        for s in split_steps:
            # base case: step is a ko, mnum, non-essential step, or '--'
            if (len(s) == 6 and s[0] == "K") or (len(s) == 6 and s[0] == "M") or (s == "--") or (
                    len(s) == 7 and s[0] == "-"):
                for p in paths_list:
                    p.extend([s])
            else:
                if s[0] == "(" and s[-1] == ")":
                    # here we try splitting to see if removing the outer parentheses will make the definition become unbalanced
                    # (the only way to figure this out is to try it because regex cannot handle nested parentheses)
                    comma_substeps = self.split_by_delim_not_within_parens(s[1:-1], ",")
                    if not comma_substeps:  # if it doesn't work, try without removing surrounding parentheses
                        comma_substeps = self.split_by_delim_not_within_parens(s, ",")
                    space_substeps = self.split_by_delim_not_within_parens(s[1:-1], " ")
                    if not space_substeps:
                        space_substeps = self.split_by_delim_not_within_parens(s, " ")
                else:
                    comma_substeps = self.split_by_delim_not_within_parens(s, ",")
                    space_substeps = self.split_by_delim_not_within_parens(s, " ")

                # complex case: no commas OR spaces outside parentheses so this is a protein complex rather than a compound step
                if len(comma_substeps) == 1 and len(space_substeps) == 1:
                    complex_components, delimiters = self.split_by_delim_not_within_parens(s, ["+", "-"],return_delims=True)
                    complex_strs = [""]

                    # reconstruct the complex (and any alternate possible complexes) while keeping the +/- structure the same
                    for i in range(len(complex_components)):
                        c = complex_components[i]
                        if c[0] == '(':
                            alts = self.split_path(c)
                            new_complex_strs = []
                            for a in alts:
                                if len(a) > 1:
                                    raise Exception
                                for cs in complex_strs:
                                    extended_complex = cs + a[0]
                                    new_complex_strs.append(extended_complex)
                            complex_strs = new_complex_strs
                        else:
                            for j in range(len(complex_strs)):
                                complex_strs[j] += c

                        if i < len(delimiters):
                            for j in range(len(complex_strs)):
                                complex_strs[j] += delimiters[i]

                    new_paths_list = []
                    for cs in complex_strs:
                        for p in paths_list:
                            p_copy = copy.copy(p)
                            p_copy.extend([cs])
                            new_paths_list.append(p_copy)
                    paths_list = new_paths_list

                # compound step case:
                else:
                    alts = self.split_path(s)
                    new_paths_list = []
                    for a in alts:
                        for p in paths_list:
                            p_copy = copy.copy(p)
                            p_copy.extend(a)
                            new_paths_list.append(p_copy)
                    paths_list = new_paths_list

        return paths_list


    def remove_non_essential_kos(self, ko_str):
        res=[]
        re_pattern=re.compile(r'-K\d{5}')
        for step in ko_str:
            temp=step
            search=re.findall(re_pattern,step)
            for s in search:
                temp=temp.replace(s,'')
            res.append(temp)
        return res

    def get_sets_module(self, string_to_search):
        module_str=string_to_search.split('hidden">')[-1].split('<br>')[0].strip().replace('<wbr>','')
        res=[]
        all_paths= self.recursive_definition_unroller(module_str)
        for i in range(len(all_paths)):
            only_essentials=self.remove_non_essential_kos(all_paths[i])
            temp=set()
            for step in only_essentials:
                step_modules=step.split('+')
                temp.update(step_modules)
            res.append(temp)
        return res


    def get_ko_from_module(self, module_id):
        url = f'https://www.genome.jp/dbget-bin/www_bget?md:{module_id}'
        print(f'Getting module {module_id}')
        webpage = None
        c = 0
        while not webpage and c <= 10:
            req = requests.get(url)
            try:
                webpage = req.text
            except Exception:
                c += 1
                continue
            print(webpage)
            start=re.search('>Definition<',webpage).span()[1]
            webpage=webpage[start:]
            end=re.search('</div></div></td></tr>',webpage).span()[0]
            webpage=webpage[:end]
            ko_str=self.get_sets_module(webpage)
        return ko_str


    def read_modules(self, modules_dict):
        tree_modules={}
        for main_path in modules_dict:
            main_path_name=main_path['name']
            if main_path_name not in tree_modules:
                tree_modules[main_path_name]={}
            sub_pathways=main_path['children']
            for sub_path in sub_pathways:
                sub_path_name=sub_path['name']
                if sub_path_name not in tree_modules[main_path_name]:
                    tree_modules[main_path_name][sub_path_name] = {}
                modules=sub_path['children']
                for module in modules:
                    module_name=module['name'].split('[')[0]
                    module_id=module_name.split()[0]
                    module_name=module_name.replace(module_id,'').strip()
                    tree_modules[main_path_name][sub_path_name][module_id]=[module_name,self.get_ko_from_module(module_id)]
        return tree_modules


    def run(self, modules_json_path: str, kegg_tree_path: str):
        modules_dict= json.load(open(modules_json_path))['children'][0]['children']
        tree_modules = self.read_modules(modules_dict=modules_dict)
        with open(kegg_tree_path, 'wb') as handle:
            pickle.dump(tree_modules, handle)


if __name__ == '__main__':
    # from https://www.genome.jp/kegg-bin/show_brite?ko00002.keg
    kegg_module = 'modules.json'
    pickle_path = 'modules.pickle'
    kegg_tree_generator = KeggTreeGenerator()
    kegg_tree_generator.run(modules_json_path=kegg_module,
                            kegg_tree_path='kegg_tree.json')

