#for consensus
from string import punctuation
from nltk import cluster
from nltk.tag.perceptron import PerceptronTagger

from nltk.corpus import wordnet as wn
from nltk.tag import SequentialBackoffTagger,_pos_tag
from nltk.probability import FreqDist
from nltk.corpus import stopwords
from nltk import edit_distance
from nltk import download as nltk_download
import pickle
import numpy
import os
import re




'''
What is the aim of this class?
To vectorize strings and measure the importance of each word within the string. This can then be used to measure similarity between strings,
How does this class work?



1- Free text pre processing where the text annexed to each annotation is split into documents:
        Pre-process strings
        Remove identifiers
        Standardize punctuation
        Remove digits that are not attached to a token
        Standardize ion patterns
        Replace Roman numerals with Arabic numerals
        Divide by parentheses 
        Unite certain terms (example: (<polymerase> <3> -> <polymerase 3>))

2- Token tagging
        pos_tag with universal token tagging (contextual)
        Wordnet token tagging (independent)
        Choose best tag (Wordnet takes priority)
        Removal of unwanted tagged tokens (determiners, pronouns, particles, and conjunctions)
        
3- Token scoring
        Try to find synonyms shared between the 2 compared documents
        Build TF-IDF vector with tokens 
        Calculate cosine distance between the two vectors
        Calculate Jaccard distance between identifiers found during pre-processing

4- If cosine distance is above the nlp_threshold consider it a match


text similarity articles:
https://towardsdatascience.com/overview-of-text-similarity-metrics-3397c4601f50
https://medium.com/@adriensieg/text-similarities-da019229c894


TO UPDATE UNIPROT LEXICON:
#https://www.uniprot.org/uniprot/?query=reviewed
go to uniprot and search for all proteins (empty search field), then choose the collumns (protein_names and annotation). Apply collumns
Select only the reviewed proteins.
Download results in tab separated format.
Delete old pickled files and set path to new file in the MANTIS.config

To update gene ontology synonyms and terms:
Go to http://geneontology.org/docs/download-ontology/ and download go.obo
Delete old pickled files and set path to new file in the MANTIS.config

'''


class Word_Weighter():

    def process_uniprot(self,line):
        if 'Entry' not in line:
            line = line.split('\t')
            header,annotation,annot_score=line
            return header,annotation

    def is_float(self,x):
        try:
            float(x)
            return True
        except:
            return False

    def is_abbreviation(self,x):
        if re.search('#[a-z]?[A-Z]+',x):
            return True
        return False


    def generate_n_grams(self,list_words):
        for n_range in self.n_grams_range:
            grams = [list_words[i:i + n_range] for i in range(len(list_words) - n_range + 1)]
            for g in grams:
                yield ' '.join(g)

    def build_frequency_dict(self):
        self.load_word_counter_pickle()
        if not self.word_counter:
            print('Building frequency dictionary with ',self.uniprot_reference,flush=True)
            stop_words = set(stopwords.words("english"))
            with open(self.uniprot_reference) as file:
                line=file.readline()
                while line:
                    line=line.strip('\n')
                    func_res=self.process_uniprot(line)
                    if func_res:
                        self.document_counter+=1
                        query,annotation=func_res
                        already_added=set()
                        #each annotation line is a document,each word a token, the whole collection of annotations is the corpus
                        list_words=self.pre_process_string(annotation)[0]
                        for sw in list_words:
                            n_grams=self.generate_n_grams(sw)
                            #for word in sw:
                            for word in n_grams:
                                if not self.is_float(word) and word not in stop_words:
                                    if word not in self.word_counter: self.word_counter[word]={'N_docs_with_token':0}
                                    if word not in already_added:
                                        self.word_counter[word]['N_docs_with_token'] += 1
                                        already_added.add(word)

                    line=file.readline()
            self.get_metrics_idfs()
            self.save_word_counter_pickle()

    def get_metrics_idfs(self):
        all_idfs=[]
        for word in self.word_counter:
            word_idf=self.calculate_idf(word)
            all_idfs.append(word_idf)
        self.min_idf=min(all_idfs)
        self.max_idf=max(all_idfs)

    def save_word_counter_pickle(self):
        with open(self.word_counter_pickled, 'wb') as handle:
            pickle.dump((self.word_counter,self.document_counter,self.min_idf,self.max_idf), handle)

    def load_word_counter_pickle(self):
        if os.path.exists(self.word_counter_pickled):
            with open(self.word_counter_pickled, 'rb') as handle:
                self.word_counter,self.document_counter,self.min_idf,self.max_idf= pickle.load(handle)

    def min_max_scale_nlp(self,X,minX,maxX):
        if minX==0 and maxX==0: return 0
        if minX==maxX: return 1
        return (X-minX)/(maxX-minX)


    #is token too common in the corpus?
    def calculate_idf(self,word):
        if word in self.words_to_remove:
            res=0
        elif word in self.word_counter:
            N_docs_with_token=self.word_counter[word]['N_docs_with_token']
            res = self.document_counter / N_docs_with_token
        elif self.is_float(word):
            res=0
        elif word in set(punctuation):
            res=0
        else:
            res=0
        return res

    #when coming from the main nlp, with all of its pre processing
    #if the word appears multiple times it should be more important?
    def calculate_tf_idf(self,list_words):
        all_words=[]
        for word in self.generate_n_grams(list_words):
            all_words.append(word)
        sentence_word_counter={}
        for word in all_words:
            if word not in sentence_word_counter: sentence_word_counter[word]=0
            sentence_word_counter[word]+=1
        res={}
        for word in sentence_word_counter:
            if word not in res: res[word]=0
            n_token_in_doc=sentence_word_counter[word]
            total_tokens_in_doc=len(all_words)
            tf=n_token_in_doc/total_tokens_in_doc
            idf = self.calculate_idf(word)
            res[word]=tf*idf
        return res

    def calculate_scaled_tf_idf(self,tf_ids):
        #we scale the weights to understand which tokens are more important within each sentence
        maxX = max(tf_ids.values())
        minX = min(tf_ids.values())
        for word in tf_ids:
            tf_ids[word] = self.min_max_scale_nlp(tf_ids[word], minX, maxX)
        return tf_ids

    def calculate_tf_idf_sentence(self, sentence):
        list_sentences_tokens = self.pre_process_string(sentence)[0]
        list_words = []
        for sentence in list_sentences_tokens:
            for word in sentence:
                list_words.append(word)
        if list_words:
            return self.calculate_tf_idf(list_words)
        else: return {}

####PREPROCESSING

class NLP_Pre_Processer():
    def remove_ecs(self, string_to_search, required_level=3):
        matches = []
        removed_ecs = set()
        # greedy match of confounders
        ec_pattern = re.compile('\d(\.(-|\d{1,3}|([a-zA-Z]\d{1,3}))){2,3}')
        search = re.finditer(ec_pattern, string_to_search)
        for i in search:
            ec = i.group()
            passed = False
            start = i.span()[0]
            end = i.span()[1]
            if len(string_to_search) > end + 1:
                if string_to_search[start - 1] != '.' and string_to_search[end - 1] != '.' and not re.match(
                        '\.|[a-zA-Z]|\d{1,3}', string_to_search[end + 1]) and not re.match('-', string_to_search[end]):
                    passed = True
            else:
                if string_to_search[start - 1] != '.':
                    passed = True
            if passed:
                if ec.count('.') >= required_level - 1:
                    if ec.count('.') + 1 - ec.count('-') >= required_level:
                        matches.append([start, end])
                        removed_ecs.add(ec)
        removed_length = 0
        current_str = str(string_to_search)
        for m in matches:
            start, end = m
            start -= removed_length
            end -= removed_length
            current_str = current_str[:start] + current_str[end:]
            removed_length += end - start
        to_remove_pattern = re.compile(
            '(\(EC\s? \)|\(EC:\s?\)|\(ec\s?\)|\(ec:\s?\)|\[EC\s?\]|\[EC:\s?\]|\[ec\s?\]|\[ec:\s?\])')
        search = re.search(to_remove_pattern, current_str)
        if search:
            current_str = current_str.replace(search.group(), '')
        return current_str, removed_ecs

    def remove_pattern(self, string_to_search, pattern):
        patterns_removed = set()
        search = re.search(pattern, string_to_search)
        while search:
            patterns_removed.add(search.group())
            start = search.span()[0]
            end = search.span()[1]
            if string_to_search[start + 1] == '(': start += 2
            if string_to_search[end - 1] == ')': end -= 1
            string_to_search = list(string_to_search)
            string_to_search[start:end] = ''
            string_to_search = ''.join(string_to_search)
            search = re.search(pattern, string_to_search)
        return string_to_search, patterns_removed

    def convert_to_arabic_digits(self, roman_digit):
        roman_numerals = [
            ('M', 1000),
            ('CM', 900),
            ('D', 500),
            ('CD', 400),
            ('C', 100),
            ('XC', 90),
            ('L', 50),
            ('XL', 40),
            ('X', 10),
            ('IX', 9),
            ('V', 5),
            ('IV', 4),
            ('I', 1)
        ]
        ix = 0
        result = 0
        while ix < len(roman_digit):
            for k, v in roman_numerals:
                if roman_digit.startswith(k, ix):
                    result += v
                    ix += len(k)
                    break
            else:
                raise ValueError('Invalid Roman number.')
        return result

    def replace_roman_numerals(self, string_to_search):
        # we wont use high roman numerals since they dont usually come up in this scenario
        roman_pattern = re.compile('[^a-zA-Z0-9][IV]+[^a-zA-Z0-9\.]')
        search = re.search(roman_pattern, string_to_search)
        while search:
            start = search.span()[0]
            end = search.span()[1]
            roman_digit = re.search('[IV]+', search.group()).group()
            string_to_search = list(string_to_search)
            converted_number = search.group().replace(roman_digit.upper(),
                                                      str(self.convert_to_arabic_digits(roman_digit)))
            string_to_search[start:end] = converted_number
            string_to_search = ''.join(string_to_search)
            search = re.search(roman_pattern, string_to_search)

        return string_to_search

    def replace_punctuation(self, sentence):
        temp_sentence = str(sentence)
        temp_sentence = temp_sentence.replace('-->', 'to')
        punctuation_set = set(punctuation)
        for i in ['\'', '-', '.', ',', '+', '(', ')', '[', ']']:
            punctuation_set.remove(i)
        punctuation_set.add(', ')
        punctuation_set.add(' - ')
        for p in punctuation_set:
            temp_sentence = temp_sentence.replace(p, ' ')
        # some terms are wrongly separated ex:GDP-4- dehydro
        temp_sentence = temp_sentence.replace('- ', '-')
        return temp_sentence

    def remove_bad_pattern_ions(self, string_to_search):
        bad_ion_pattern = re.compile('\(\d\s\)')
        search = re.search(bad_ion_pattern, string_to_search)
        while search:
            ion_str = search.group()
            new_ion_str = ion_str.replace(' ', '')
            string_to_search = string_to_search.replace(ion_str, new_ion_str)
            search = re.search(bad_ion_pattern, string_to_search)
        return string_to_search

    def simplify_ions(self, string_to_search):
        ion_pattern = re.compile('(\(\d\([\+\-]\)\))|(\(\d\)\([\+\-])\)')
        search = re.search(ion_pattern, string_to_search)
        while search:
            ion_str = search.group()
            new_ion_str = '(' + ion_str.replace('(', '').replace(')', '') + ')'
            string_to_search = string_to_search.replace(ion_str, new_ion_str)
            search = re.search(ion_pattern, string_to_search)
        return string_to_search

    def get_parentheses_above_range(self, starting_range, p2_ranges):
        res = []
        for p2 in p2_ranges:
            if p2[0] > starting_range[0]: res.append(p2)
        return res

    def get_parentheses_pairs(self, p1_ranges, p2_ranges):
        res = []
        remaining = []

        if not p1_ranges or not p2_ranges: return res, remaining
        while p1_ranges:
            current_p1 = p1_ranges.pop(-1)
            p2_above_range = self.get_parentheses_above_range(current_p1, p2_ranges)
            if not p2_above_range:
                remaining = list(p1_ranges)
                break
            p2_for_p1 = p2_above_range[0]
            p2_ranges.remove(p2_for_p1)
            res.append([current_p1, p2_for_p1])
        return res, remaining

    def remove_parentheses(self, string_to_search):
        p1_pattern = re.compile('\(')
        p2_pattern = re.compile('\)')
        p1_search = re.finditer(p1_pattern, string_to_search)
        p2_search = re.finditer(p2_pattern, string_to_search)
        p1_search = [i.span() for i in p1_search]
        p2_search = [i.span() for i in p2_search]
        res = list(string_to_search)
        # whe parentheses are loose we just removed them
        if not p2_search:
            for p1 in p1_search:
                res[p1[0]] = ''
        if not p1_search:
            for p2 in p2_search:
                res[p2[0]] = ''
        p_pairs, remaining = self.get_parentheses_pairs(p1_search, p2_search)
        for r in remaining:
            res[r[0]] = ''
        to_remove = []
        for i in p_pairs:
            if string_to_search[i[0][0] - 1:i[0][1]] == ' (' and (string_to_search[i[1][0]:i[1][1] + 1] == ') '
                                                                  or string_to_search[i[1][0]:i[1][1] + 2] == '). '
                                                                  or string_to_search[i[1][0]:i[1][1] + 2] == ').\t'
                                                                  or string_to_search[i[1][0]:i[1][1] + 1] == ')\t'):
                to_remove.append(i)
        for t in to_remove:
            res[t[0][0]] = '#SPLIT#'
            res[t[1][0]] = '#SPLIT#'
        res = ''.join(res)
        return res

    def unite_terms(self, string_to_process):
        string_list = string_to_process.split()
        res = []
        c = 0
        for i in range(len(string_list)):
            if i == 0:
                res.append(string_list[i])
            else:
                if (len(string_list[i]) == 1 and string_list[i].isupper()) or re.search('\d+\s', string_list[i]):
                    c += 1
                    res[i - c] += ' ' + string_list[i]
                else:
                    res.append(string_list[i])
        c = 0
        res=[word for word in res if word]
        for i in range(len(res)):
            if res[i][0] == '-':
                res[i] = res[i].lstrip('-')
                c += 1
            if res[i][-1] == '-':
                res[i] = res[i].rstrip('-')
                c += 1
            if (res[i].count('(') == 1 and res[i].count(')') == 0) or \
                    (res[i].count(')') == 1 and res[i].count('(') == 0):
                c += 1
                res[i] = res[i].strip(')')
        return res

    def get_token_to_merge(self, list_of_tokens, to_add_pos):
        for i in range(len(list_of_tokens)):
            if i > to_add_pos:
                passed = True
                for t in list_of_tokens[::-1][i]:
                    if '#' in t: passed = False
                if passed:
                    list_of_tokens[::-1][i].extend(list_of_tokens[::-1][to_add_pos][1:])
                    to_remove = list_of_tokens[::-1].pop(to_add_pos)
                    list_of_tokens.remove(to_remove)
                    return list_of_tokens

    def connect_gapped_token(self, list_of_tokens):
        to_append = []
        for i in range(len(list_of_tokens)):
            if list_of_tokens[i][0] == '#':
                to_append.append(i)
        for i in to_append[::-1]:
            self.get_token_to_merge(list_of_tokens, len(list_of_tokens) - i - 1)
        for li in range(len(list_of_tokens)):
            for ti in range(len(list_of_tokens[li])):
                if not self.is_abbreviation(list_of_tokens[li][ti]):
                    list_of_tokens[li][ti] = list_of_tokens[li][ti].replace('#.', '')
                    list_of_tokens[li][ti] = list_of_tokens[li][ti].replace('#', '')
                else:
                    list_of_tokens[li][ti]=''
        res = []
        for lt in list_of_tokens:
            temp = []
            for token in lt:
                if token:
                    temp.append(token)
            if temp:
                res.append(temp)
        return res

    def final_processing(self, string_to_process):
        new_str = string_to_process.strip(' ')
        new_str=new_str.replace('[]','')
        new_str=new_str.replace('()','')
        new_str = new_str.replace('. ', '!NEWLINE!')
        new_str = new_str.replace('\t', '!NEWLINE!')
        lines = new_str.split('!NEWLINE!')
        res = []

        for line in lines:
            split_parentheses = line.split('#SPLIT')
            for current_str in split_parentheses:
                current_str = self.unite_terms(current_str)
                #current_str = [i.lower() for i in current_str if not self.is_float(i)]
                if current_str:
                    res.append(current_str)
            if 'SPLIT' in string_to_process:
                res = self.connect_gapped_token(res)
        return res

    def remove_go_obo_identifiers(self, sentence):
        '''
        ids to keep
        vz =  https://viralzone.expasy.org/
        SO = http://www.sequenceontology.org/
        metacyc
        reactome
        hgnc
        pfam
        chebi
        brenda

        cant keep KEGG because GO.obo has no distinction between KEGG's identifiers types...
        '''
        go_obo_pattern = re.compile('\[('
                                    'goc|PR|CL|Wikipedia|CORUM|MetaCyc|GOC|ISBN|PMID|Reactome|CHEBI|GO|VZ|vz|gOC|HGNC|KEGG|KEGG_REACTION|UBERON|Pfam|RESID|MA|SO|UniPathway|MP|BRENDA|DOI|pmid|Wikipeda|Wikilpedia|MGI|DDANAT|PO|ABA'
                                    '):[A-Za-z\d\-]+(,\s('
                                    'goc|PR|CL|Wikipedia|CORUM|MetaCyc|GOC|ISBN|PMID|Reactome|CHEBI|GO|VZ|vz|gOC|HGNC|KEGG|KEGG_REACTION|UBERON|Pfam|RESID|MA|SO|UniPathway|MP|BRENDA|DOI|pmid|Wikipeda|Wikilpedia|MGI|DDANAT|PO|ABA'
                                    '):[A-Za-z\d\-]+)*\]')
        res, go_obo_ids = self.remove_pattern(sentence, go_obo_pattern)
        go_obo_ids_res = {}
        http_pattern = re.compile('\[http.*\]')
        search_http = re.search(http_pattern, res)
        if search_http:
            res = res.replace(search_http.group(), '')
        if go_obo_ids:
            go_obo_ids = go_obo_ids.pop()
            go_obo_ids = go_obo_ids.replace('[', '')
            go_obo_ids = go_obo_ids.replace(']', '')
            go_obo_ids = go_obo_ids.split(', ')
            go_obo_ids_res = {}
            for i in go_obo_ids:
                go_obo_type, go_obo_id = i.split(':')
                if go_obo_type in ['vz', 'VZ', 'SO', 'MetaCyc', 'Reactome', 'HGNC', 'Pfam', 'CHEBI', 'BRENDA']:
                    if go_obo_type == 'vz':
                        go_obo_type = 'viralzone'
                    elif go_obo_type == 'VZ':
                        go_obo_type = 'viralzone'
                    elif go_obo_type == 'SO':
                        go_obo_type = 'seq_onto'
                    elif go_obo_type == 'BRENDA':
                        go_obo_type = 'enzyme_ec'
                    else:
                        go_obo_type = go_obo_type.lower()
                    if go_obo_type not in go_obo_ids_res:
                        go_obo_ids_res[go_obo_type] = set()
                    go_obo_ids_res[go_obo_type].add(go_obo_id)
        return res, go_obo_ids_res

    def remove_common_identifiers(self, sentence):
        res = ' ' + str(sentence) + ' '
        dict_ids = {'enzyme_ec': set(),
                    'tcdb': set(),
                    'kegg_ko': set(),
                    'tigrfam': set(),
                    'pfam': set(),
                    'cog': set(),
                    'go': set(),
                    'others': set()
                    }
        res, removed_ecs = self.remove_ecs(res, required_level=1)
        res = res.replace(' enzyme_ec:', '')
        dict_ids['enzyme_ec'].update(removed_ecs)
        tcdb_pattern = re.compile('\(TC\s\d\.[A-Z\-](\.(\d+|\-)){1,2}\)')
        ko_pattern = re.compile('K\d{4,}')
        tigrfam_pattern = re.compile('TIGR\d+')
        duf_pattern = re.compile('(DUF|duf)\d+')
        pfam_pattern = re.compile('((U|u)?PF|pf)\d+')
        cog_pattern = re.compile('(COG|cog)\d+')
        go_pattern = re.compile('GO:?\d+')
        res, removed_ids = self.remove_pattern(res, tcdb_pattern)
        res = res.replace(' tcdb:', '')
        dict_ids['tcdb'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, tigrfam_pattern)
        res = res.replace(' tigrfam:', '')
        dict_ids['tigrfam'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, ko_pattern)
        res = res.replace(' kegg_ko:', '')
        dict_ids['kegg_ko'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, duf_pattern)
        res = res.replace(' pfam:', '')
        dict_ids['pfam'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, pfam_pattern)
        dict_ids['pfam'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, cog_pattern)
        res = res.replace(' cog:', '')
        dict_ids['cog'].update(removed_ids)
        res, removed_ids = self.remove_pattern(res, go_pattern)
        res = res.replace(' go:', '')
        removed_ids = {re.search('\d+', i).group() for i in removed_ids}
        dict_ids['go'].update(removed_ids)
        res, go_obo_identifiers = self.remove_go_obo_identifiers(res)
        for go_obo_type in go_obo_identifiers:
            if go_obo_type == 'enzyme_ec':
                dict_ids['enzyme_ec'].update(go_obo_identifiers[go_obo_type])
            else:
                dict_ids[go_obo_type] = go_obo_identifiers[go_obo_type]
        return res, dict_ids

    def pre_process_string(self, sentence):
        if not sentence: return [], []
        res,dict_ids = self.remove_common_identifiers(sentence)
        res = self.replace_punctuation(res)
        digit_pattern = re.compile('[^A-Za-z][\s\(](\d+(\.\d+)?)[\s\)]')
        res, _ = self.remove_pattern(res, digit_pattern)

        id_pattern = re.compile('[A-Z]+\d{3,}')
        res, ids_removed = self.remove_pattern(res, id_pattern)
        dict_ids['others'].update(ids_removed)
        # cleaning blank parenthesis
        res, _ = self.remove_pattern(res, '\(\s+\)')
        res, _ = self.remove_pattern(res, '\[\s+\]')
        res = self.remove_bad_pattern_ions(res)
        res = self.simplify_ions(res)
        res = self.replace_roman_numerals(res)
        res = self.remove_parentheses(res)
        all_ids = set()
        for id_type in dict_ids:
            if id_type == 'others':
                all_ids.update(dict_ids['others'])
            else:
                for id_str in dict_ids[id_type]:
                    all_ids.add(id_type + ':' + id_str)
        return self.final_processing(res), all_ids



class WordNetTagger(SequentialBackoffTagger):
    def __init__(self,perceptron_tagger,go_terms=None, *args, **kwargs):
        SequentialBackoffTagger.__init__(self, *args, **kwargs)
         #for universal tagger
        self.placeholder = 'XXXXX'
        self.perceptron_tagger=perceptron_tagger
        self.go_terms=set()
        for g in go_terms:
            #these tags were manually reviewed so they wouldnt intefere with go terms tagging
            if self.tag_tokens_perceptron([g]) not in ['ADP','CONJ','DET','PRON','PRT']:
                self.go_terms.add(g)
        self.wordnet_tag_map = {'NOUN': 'NOUN','ADJ': 'ADJ',
                                'ADV': 'ADV','VERB': 'VERB',
                                #placeholder for entities
                                self.placeholder:self.placeholder}

    def tag_tokens_perceptron(self, tokens):
        return _pos_tag(tokens, tagset='universal', tagger=self.perceptron_tagger,lang='eng')

    def most_frequent(self,list_to_test):
        return max(set(list_to_test), key=list_to_test.count)

    def choose_tag(self, tokens, index, history):
        word = tokens[index]
        word=word.strip()
        if word ==self.placeholder: return self.placeholder
        #adding terms from pre-processing
        if len(word.split(' '))>1:
            terms=word.split(' ')
            word=terms[0]
        fd = FreqDist()
        for synset in wn.synsets(word):
            lex_name=synset.lexname().split('.')[0]
            fd[lex_name] += 1
        wordnet_res,go_res=None,None
        if fd.keys(): wordnet_res= self.wordnet_tag_map.get(fd.max())
        tagger_res=self.tag_tokens_perceptron([word])[0][1]
        if word in self.go_terms: go_res= 'NOUN'
        if wordnet_res: return wordnet_res
        elif tagger_res:
            if tagger_res!='NOUN':
                return tagger_res
        elif go_res: return go_res
        else: return None

class MANTIS_NLP(NLP_Pre_Processer,Word_Weighter):
    def __init__(self,n_grams_range=[1]):
        self.download_nltk_resources()
        self.tagger = PerceptronTagger()
        self.n_grams_range=n_grams_range
        self.set_path_go_terms_nlp()
        self.set_path_uniprot_proteins_nlp()
        str_n_gram = '_'.join([str(i) for i in self.n_grams_range])
        self.word_counter_pickled = self.uniprot_reference+'_n_grams_'+str_n_gram+'.pickle'
        self.word_counter={}
        self.document_counter=0

        self.good_identifiers={'enzyme_ec','tcdb','kegg_ko','tigrfam','pfam','cog','go','viralzone','seq_onto'}

        self.words_to_remove = ['mainrole','sub1role','protein','proteins',
                                'enzyme','enzymes','putative','activity',
                                'process','unknown','function','functions',
                                'processes'
                                ]
        self.build_frequency_dict()
        self.pickled_go_syns = self.go_terms_path+'.pickle_syns'
        self.pickled_go_terms = self.go_terms_path+'.pickle_terms'
        self.pickled_go_dict = self.go_terms_path+'.pickle_dict'
        self.go_syns=set()
        self.go_terms=set()
        self.go_dict=dict()
        self.parse_go_terms()
        self.wordnet_tagger=WordNetTagger(go_terms=self.go_terms,perceptron_tagger=self.tagger)
        self.tags={}
        self.nlp_threshold=0.8

    def tag_tokens_perceptron(self, tokens):
        return _pos_tag(tokens, tagset='universal', tagger=self.tagger,lang='eng')

    def tag_tokens_wordnet(self, tokens):
        return  self.wordnet_tagger.tag(tokens)



    def __str__(self):
        res='Tags:\n'
        for t in self.tags:
            res+=t+': '+str(self.tags[t])+'\n'
        return res

    def download_nltk_resources(self):
        try:
            nltk_download('stopwords',quiet=True)
        except:
            print('Already downloaded stopwords')
        try:
            nltk_download('averaged_perceptron_tagger',quiet=True)
        except:
            print('Already downloaded Perceptron tagger')
        try:
            nltk_download('universal_tagset',quiet=True)
        except:
            print('Already downloaded Universal tagset!')
        try:
            nltk_download('wordnet',quiet=True)
        except:
            print('Already downloaded Wordnet!')

    ####GO SCORING
    def save_go_pickle(self):
        with open(self.pickled_go_syns, 'wb') as handle:
            pickle.dump(self.go_syns, handle,protocol=4)
        with open(self.pickled_go_terms, 'wb') as handle:
            pickle.dump(self.go_terms, handle,protocol=4)
        with open(self.pickled_go_dict, 'wb') as handle:
            pickle.dump(self.go_dict, handle,protocol=4)

    def load_go_pickle(self):
        if os.path.exists(self.pickled_go_syns):
            with open(self.pickled_go_syns, 'rb') as handle:
                self.go_syns = pickle.load(handle)
        if os.path.exists(self.pickled_go_terms):
            with open(self.pickled_go_terms, 'rb') as handle:
                self.go_terms = pickle.load(handle)
        if os.path.exists(self.pickled_go_dict):
            with open(self.pickled_go_dict, 'rb') as handle:
                self.go_dict = pickle.load(handle)

    def token_list_too_small(self,tokens_list,all_tokens_list):
        res=[]
        if len(tokens_list)>1: return False
        for tl in all_tokens_list:
            res.append(len(tl))
        max_len= max(res)
        if len(tokens_list)==1 and len(tokens_list)<max_len:
            return True
        return False

    def parse_go_terms(self):
        self.load_go_pickle()
        if not self.go_syns or not self.go_terms or not self.go_dict:
            print('Parsing GO terms with ',self.go_terms_path,flush=True)
            go_terms=set()
            go_dict={}
            with open(self.go_terms_path) as file:
                line=file.readline()
                while line:
                    line=line.strip('\n')
                    if line:
                        if 'id: ' in line[0:4]:
                            go_id=line.replace('id: ','')
                            go_id=go_id.lower()
                            go_dict[go_id]={'synonyms':set(),'identifiers':set()}
                        elif 'alt_id: ' in line[0:8]:
                            alt_go_id=line.replace('alt_id: ','')
                            alt_go_id=alt_go_id.lower()
                            go_dict[alt_go_id]=go_dict[go_id]


                        elif 'name: ' in line[0:6]:
                            go_name=line.replace('name: ','')
                            split_tokens,all_ids=self.pre_process_string(go_name)
                            for i in split_tokens:
                                tags=self.tag_tokens_perceptron(i)
                                for t in tags:
                                    if t[1] not in ['ADP','CONJ','DET','PRON','PRT'] and\
                                        t[0] not in self.words_to_remove and\
                                        len(t[0])>1 and\
                                        not  re.search('\d+',t[0]):
                                        go_terms.add(t[0])
                            for token_list in split_tokens:
                                if not self.token_list_too_small(token_list,split_tokens):
                                    go_dict[go_id]['synonyms'].add(' '.join(token_list))
                            go_dict[go_id]['identifiers'].update(all_ids)
                        elif 'synonym: ' in line[0:9] and 'EXACT' in line:
                            go_syn=line.replace('synonym: ','')
                            go_syn=go_syn.replace('EXACT ','')
                            go_syn=go_syn.replace('[]','')

                            split_tokens,all_ids=self.pre_process_string(go_syn)
                            for token_list in split_tokens:
                                if not self.token_list_too_small(token_list,split_tokens):
                                    go_dict[go_id]['synonyms'].add(' '.join(token_list))
                            go_dict[go_id]['identifiers'].update(all_ids)

                    line = file.readline()
                    if '[Typedef]' in line: line=None
            for go_id in go_dict:
                self.go_syns.add(frozenset(go_dict[go_id]['synonyms']))
            self.go_terms=go_terms
            self.go_dict=go_dict
        self.save_go_pickle()

    def has_go_match(self,test_syn,ref_syn):
        for syn_set in self.go_syns:
            if ref_syn in syn_set and test_syn in syn_set:
                return True
        return False



    ####NLP SCORING

    def remove_unwanted_tokens(self,tagged_tokens):
        '''
        POS tag list:
        for the universal tagset:
            ADJ 	adjective 	new, good, high, special, big, local
            ADP 	adposition 	on, of, at, with, by, into, under
            ADV 	adverb 	really, already, still, early, now
            CONJ 	conjunction 	and, or, but, if, while, although
            DET 	determiner, article 	the, a, some, most, every, no, which
            NOUN 	noun 	year, home, costs, time, Africa
            NUM 	numeral 	twenty-four, fourth, 1991, 14:24
            PRT 	particle 	at, on, out, over per, that, up, with
            PRON 	pronoun 	he, their, her, its, my, I, us
            VERB 	verb 	is, say, told, given, playing, would
            . 	punctuation marks 	. , ; !
            X 	other 	ersatz, esprit, dunno, gr8, univeristy
        '''
        res = []
        #ideally we'd only keep the nouns and numbers (maybe verbs?)but there's a lot of false positives...
        tags_to_remove=[self.wordnet_tagger.placeholder,'DET','PRON','PRT','CONJ']
        #from tigrfam
        stop_words = set(stopwords.words("english"))
        for t in tagged_tokens:
            if t[1] not in ['NOUN']:
                if t[1] not in self.tags: self.tags[t[1]]=set()
                self.tags[t[1]].add(t[0])
            if t[1] not in tags_to_remove and len(t[0])>1 and t[0] not in self.words_to_remove and t[0] not in stop_words:
                res.append(t[0])
        return res

    def choose_best_tagging(self,wornet_tagging,default_tagging):
        res=[]
        for word_tag_i in range(len(wornet_tagging)):
            if default_tagging[word_tag_i][1]=='NUM':res.append(default_tagging[word_tag_i])
            elif wornet_tagging[word_tag_i][1]: res.append(wornet_tagging[word_tag_i])
            else: res.append(default_tagging[word_tag_i])
        return res

    def process_string_nlp(self,tokens_lists):
        res=[]
        for tokens in tokens_lists:
            # word lemmatization - doesnt make much sense since we are not trying to classify text meaning. It won't change the tokens that much and might even produce errors
            # word stemming - makes a bit more sense but once again it also produces weird errors
            temp_tokens=[i for i in tokens if i]
            default_tagged_tokens = self.tag_tokens_perceptron(temp_tokens)
            _ = self.tag_tokens_wordnet(temp_tokens)
            wordnet_tagged_tokens = self.wordnet_tagger.tag(temp_tokens)
            tagged_tokens=self.choose_best_tagging(wordnet_tagged_tokens,default_tagged_tokens)
            removed_tags = self.remove_unwanted_tokens(tagged_tokens)
            res.append(removed_tags)
        return res


    def get_best_syn(self,original_syns_set,original_list_words,edit_dist_perc=0.05):
        list_words=set([i.lower() for i in original_list_words])
        syns_set=set([i.lower() for i in original_syns_set])
        for syn in syns_set:
            if syn in list_words: return syn
        for syn in syns_set:
            for word in list_words:
                max_len=max([len(syn),len(word)])
                if edit_distance(syn,word)/max_len<=edit_dist_perc:
                    return word


    #maybe wont make much of a difference since the most significant words wont be in the dictionary anyway
    def improve_syn_match(self,vector1,vector2):
        #improving match by replacing vector1 words by synonyms that exist in vector2
        res=[]
        for w in vector1:
            syn_set={w}
            for synset in wn.synsets(w):
                for lemma in synset.lemmas():
                    syn_set.add(lemma.name())
            best_syn=self.get_best_syn(syn_set,vector2)
            if best_syn: res.append(best_syn)
            else: res.append(w)
        return res

    def remove_plurals(self,vector1,vector2):
        res=[]
        for w in vector1:
            syn_set={w}
            w2=None
            #latin plurals
            if w.endswith('ae'): w2=w[:-1]
            elif w.endswith('exes'): w2=w[:-2]
            elif w.endswith('ices'): w2=w[:-4]+'ex'
            elif w.endswith('eaus'): w2=w[:-1]
            elif w.endswith('eaux'): w2=w[:-1]
            elif w.endswith('ia'): w2=w[:-1]+'on'
            elif w.endswith('ions'): w2=w[:-1]
            elif w.endswith('es'): w2=w[:-2]+'is'
            elif w.endswith('os'): w2=w[:-1]
            elif w.endswith('oes'): w2=w[:-2]
            elif w.endswith('uses'): w2=w[:-2]
            elif w.endswith('i'): w2=w[:-1]+'us'
            #normal plurals
            elif w.endswith('s'): w2=w[:-1]
            if w2: syn_set.add(w2)
            #other latin plurals
            if w.endswith('a'):
                w3=w[:-1]+'um'
                syn_set.add(w3)
            if w.endswith('a'):
                w3=w[:-1]+'on'
                syn_set.add(w3)
            if w.endswith('ic'):
                w3=w[:-1]+'on'
                syn_set.add(w3)
            if w.endswith('i'):
                w3=w[:-1]+'on'
                syn_set.add(w3)
            best_syn=self.get_best_syn(syn_set,vector2)
            if best_syn: res.append(best_syn)
            else: res.append(w)
        return res

    def build_vector(self,iterable1, iterable2):
        counter1 = self.calculate_tf_idf(iterable1)
        counter1= self.calculate_scaled_tf_idf(counter1)
        counter2 = self.calculate_tf_idf(iterable2)
        counter2= self.calculate_scaled_tf_idf(counter2)
        #print('counters',counter1,counter2)
        all_items = set(counter1.keys()).union(set(counter2.keys()))
        vector1,vector2=[],[]
        for word in all_items:
            if word in counter1:vector1.append(counter1[word])
            else: vector1.append(0)
            if word in counter2:vector2.append(counter2[word])
            else: vector2.append(0)
        return vector1, vector2

    def jaccard_distance(self,label1, label2):
        union_labels=len(label1.union(label2))
        intersection_labels=len(label1.intersection(label2))
        if not union_labels: return 1
        return (union_labels - intersection_labels) / union_labels

    def score_text(self,vector_1,vector_2):
        if not vector_1 or not vector_2: return 0
        if self.wordnet_tagger.placeholder.lower() in vector_1: vector_1.remove(self.wordnet_tagger.placeholder.lower())
        if self.wordnet_tagger.placeholder.lower() in vector_2: vector_2.remove(self.wordnet_tagger.placeholder.lower())
        #improve_syn_match adds tokens in lowercase, so we need to convert vector 2 to lower case as well (also the corpus is in lowercase so thats the only way to get tf-idf)
        improved_vector_1 = self.remove_plurals(vector_1,vector_2)
        #improved_vector_1 = [t.lower() for t in vector_1]
        improved_vector_2 = [t.lower() for t in vector_2]
        v1,v2=self.build_vector(improved_vector_1,improved_vector_2)
        if not any(v1) or not any(v2):            text_score=0
        else:            text_score=1-cluster.util.cosine_distance(v1,v2)
        return text_score

    def check_token_intersection(self,tokens_lists_1,tokens_lists_2):
        for t1 in tokens_lists_1:
            for t2 in tokens_lists_2:
                if set(t1).intersection(set(t2)): return True
        return False


    def get_similarity_score(self,string_1, string_2,stdout_file=None):
        if not string_1 or not string_2: return -1
        temp_string_1=str(string_1)
        temp_string_2=str(string_2)
        #NLP SCORING
        tokens_lists_1,ids_removed_1 = self.pre_process_string(temp_string_1)
        tokens_lists_2,ids_removed_2 = self.pre_process_string(temp_string_2)
        token_intersection = self.check_token_intersection(tokens_lists_1,tokens_lists_2)
        if not token_intersection: return 0.0
        vector_1_list = self.process_string_nlp(tokens_lists_1)
        vector_2_list = self.process_string_nlp(tokens_lists_2)
        text_score=0

        for vector_1 in vector_1_list:
            for vector_2 in vector_2_list:
                current_score=self.score_text(vector_1,vector_2)
                text_score+=current_score
                #if current_score:
                    #print('#####',vector_1,'#####',vector_2,'#####',current_score,'#####', flush=True, file=stdout_file)
        good_identifiers_1 = {i for i in ids_removed_1 if i.split(':')[0] in self.good_identifiers}
        good_identifiers_2 = {i for i in ids_removed_2 if i.split(':')[0] in self.good_identifiers}
        if good_identifiers_1.intersection(good_identifiers_2): return 1

        ids_removed_score= 1-self.jaccard_distance(set(ids_removed_1),set(ids_removed_2))
        #print('text_score',text_score,vector_1,vector_2)
        #print('string_1',temp_string_1)
        #print('string_2',temp_string_2)
        #more weight for the entities
        if ids_removed_score:
            score=(2*ids_removed_score+text_score)/3
        else: score=text_score
        #GO SCORING
        if score<self.nlp_threshold and score >0.3:
            if self.has_go_match(string_1.lower(),string_2.lower()): go_score=1
            else: go_score=0
            score+=go_score
            score/=2
        return score





if __name__ == '__main__':
    #nlp=MANTIS_NLP()
    #nlp.go_terms_path='/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/Resources/Gene_Ontology/go.obo'
    #nlp.uniprot_reference='/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/Resources/Uniprot/uniprot-filtered-reviewed_yes.tab'
    s='biotin:apo-acetyl-CoA:carbon-dioxide ligase (ADP-forming) ligase (AMP-forming)'
    #nlp.false_init()
    #print(nlp.pre_process_string(s))
    #nlp.parse_go_terms()
    tool='Citrase (pro-3S)-lyase alpha chain'
    ref = 'Citrase alpha chain'
    #print(nlp.calculate_tf_idf_sentence('Unknown function (mainrole)'))
    here=['ATP-binding', 'protein', 'IstB', 'SW', 'ISTB', 'ECOLI', 'aa', 'fasta', 'scores', 'E(47.4', 'id', 'in', '249', 'aa']

    from nltk.tag import pos_tag

    a=pos_tag(here, tagset='universal')
    print(a)