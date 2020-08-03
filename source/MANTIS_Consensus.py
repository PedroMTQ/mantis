try:
    from source.cython_src.get_non_overlapping_hits import get_non_overlapping_hits
    from source.MANTIS_Assembler import *
except:
    from cython_src.get_non_overlapping_hits import get_non_overlapping_hits
    from MANTIS_Assembler import *

'''
This class creates a consensus for the integrated annotations. 

It uses NLP to vectorize and score the free text in the annotations and compares IDs between data sources, if above a certain threshold, the different sources are considered to be in consensus

After finding a consensus it also scores the hits consensus (ML) depending on:
-hit-evalue (average per consensus source)
-number of consensus annotations - non consensual annotations
-summed scored annotations/number sources
-existence of identifiers among consensus annotations


'''

class MANTIS_Consensus(MANTIS_NLP):

    def __init__(self):
        MANTIS_NLP.__init__(self)


    def get_hmm_weight(self,hmm):
        if hmm in self.mantis_hmm_weights:
            return self.mantis_hmm_weights[hmm]
        else: return self.mantis_hmm_weights['else']

    def add_from_go_obo(self,query_dict):
        for hmm_file in query_dict:
            for hmm_hit in query_dict[hmm_file]:
                temp_hmm_hit_dict_identifiers=set(query_dict[hmm_file][hmm_hit]['identifiers'])
                temp_hmm_hit_dict_description=set(query_dict[hmm_file][hmm_hit]['description'])
                for id_str in query_dict[hmm_file][hmm_hit]['identifiers']:
                    if 'go:' in id_str[0:3]:
                        if id_str in self.go_dict:
                            extra_info=self.go_dict[id_str]
                            for go_id in extra_info['identifiers']:
                                temp_hmm_hit_dict_identifiers.add(go_id)
                            for syn in extra_info['synonyms']:
                                temp_hmm_hit_dict_description.add(syn)
                #we use the temp info to help matching annotations, but we dont want to add them to the final annotations
                query_dict[hmm_file][hmm_hit]['temp_identifiers']=set(temp_hmm_hit_dict_identifiers)
                for id_str in temp_hmm_hit_dict_identifiers:
                    if 'enzyme_ec:' in id_str and '-' in id_str:
                        query_dict[hmm_file][hmm_hit]['temp_identifiers'].remove(id_str)
                query_dict[hmm_file][hmm_hit]['temp_description']=temp_hmm_hit_dict_description

    def read_interpreted_annotation(self,interpreted_annotation_tsv):
        with open(interpreted_annotation_tsv) as file:
            line=file.readline()
            line=file.readline()
            dict_annotations={}
            while line:
                line=line.strip('\n')
                line=line.strip()
                line=line.split('\t')
                query, hmm_file,hmm_hit,evalue,query_len,query_start,query_end =line[0],line[1],line[2],line[4],line[5],line[6],line[7]
                annotations=line[11:]
                if query not in dict_annotations: dict_annotations[query]={'query_len':int(query_len),'hmm_files':{}}
                #since the same hmm file can have different hits for the same query
                if hmm_file not in dict_annotations[query]['hmm_files']:
                    dict_annotations[query]['hmm_files'][hmm_file]= {}
                dict_annotations[query]['hmm_files'][hmm_file][hmm_hit]={'identifiers': set(), 'description': set(),
                                                            'evalue': float(evalue),
                                                            'query_start': int(query_start), 'query_end': int(query_end)}
                for annot in annotations:
                    if 'description:' in annot[0:len('description:')]:
                        annot_text=annot.split(':')[1:]
                        annot_text=' '.join(annot_text)
                        dict_annotations[query]['hmm_files'][hmm_file][hmm_hit]['description'].add(annot_text)
                    elif 'is_essential_gene:' in annot:
                        dict_annotations[query]['is_essential_gene']=True
                    else:
                        dict_annotations[query]['hmm_files'][hmm_file][hmm_hit]['identifiers'].add(annot)
                line=file.readline()
        for query in dict_annotations:
           self.add_from_go_obo(dict_annotations[query]['hmm_files'])
        return dict_annotations


    def generate_consensus(self,interpreted_annotation_tsv,stdout_file_path=None):
        dict_annotations=self.read_interpreted_annotation(interpreted_annotation_tsv)
        self.query_match_hits = {}
        for query in dict_annotations:
            hits = []
            query_dict=dict_annotations[query]
            query_len = query_dict['query_len']
            hmm_files = query_dict['hmm_files']
            for hmm_file in hmm_files:
                for hmm_hit in hmm_files[hmm_file]:
                    hits.append([hmm_file, hmm_hit, hmm_files[hmm_file][hmm_hit]])
            if self.domain_algorithm=='heuristic':
                best_hits = self.get_best_hits_approximation_Consensus(list(hits))
            elif self.domain_algorithm=='lowest':
                best_hits = self.get_lowest_hit_Consensus(list(hits))
            else:
                best_hits=self.get_best_cython_consensus(list(hits),query_len,query,stdout_file_path=stdout_file_path)
            consensus_hits,total_hits,hmm_files_consensus,hmm_names_consensus=self.add_to_best_hits(best_hits,dict_annotations[query]['hmm_files'])
            if 'is_essential_gene' in dict_annotations[query]:  is_essential=True
            else:                                               is_essential=False
            consensus_line=self.generate_consensus_line(query,dict_annotations[query]['hmm_files'],is_essential,consensus_hits,total_hits,hmm_files_consensus,hmm_names_consensus)
            yield consensus_line


    #same method as non_overlapping hits from MANTIS_Processor , with some minor alterations to account for consensus
    #@timeit_function
    def get_best_cython_consensus(self,hits,query_len,query,stdout_file_path=None):
        try:
            best_hits = self.get_best_hits_Consensus(list(hits),query_len,time_limit=60)
        except (TimeoutError,RecursionError):
            best_hits = self.get_best_hits_approximation_Consensus(list(hits))
            stdout_file=open(stdout_file_path,'a+')
            print('Query '+query+' was approximated during consensus generation', flush=True, file=stdout_file)
            stdout_file.close()
        return best_hits

    def get_lowest_hit_Consensus(self,query_hits):
        '''
        this will take the hit with the lowest evalue, regardless if there are multiple domains
        '''
        if not query_hits: return []
        query_hits = self.sort_by_evalue_Consensus(query_hits)
        lowest_hit= query_hits.pop(0)
        combo=[lowest_hit]
        return combo


    def query_hits_to_cython_Consensus(self, query_hits):
        conversion_dict = {}
        res = set()
        for hit_i in range(len(query_hits)):
            hmm_file=query_hits[hit_i][0]
            hmm_hit=query_hits[hit_i][1]
            hit_info=query_hits[hit_i][2]
            query_start=hit_info['query_start']
            query_end=hit_info['query_end']
            res.add(tuple([
                hit_i,
                query_start,
                query_end,
                hmm_hit
            ]))
            conversion_dict[hit_i] = [hmm_file,hmm_hit,hit_info]
        return res, conversion_dict

    def cython_to_query_hits_Consensus(self, cython_hits, conversion_dict):
        res = []
        for hit in cython_hits:
            res.append(conversion_dict[hit[0]])
        return res

    def get_min_max_Consensus(self,cython_possible_combos,conversion_dict):
        #highest evalue
        minX_evalue=-999
        #lowest evalue
        maxX_evalue=999
        for len_combo in cython_possible_combos:
            for cython_combo in cython_possible_combos[len_combo]:
                combo=self.cython_to_query_hits_Consensus(cython_combo,conversion_dict)
                for hit in combo:
                    hmm_file,hmm_hit,hit_info=hit
                    #we dont consider 0 for the scaling, 0 will always be scaled to max/1
                    if hit_info['evalue']:
                        evalue = log10(hit_info['evalue'])
                        if evalue>minX_evalue: minX_evalue=evalue
                        if evalue<maxX_evalue:maxX_evalue=evalue
        return minX_evalue,maxX_evalue



    def get_annotation_score(self,hit_info):
        if hit_info['temp_identifiers'] and hit_info['temp_description']:
            annotation_score = 1
        elif hit_info['temp_identifiers']:
            annotation_score = 0.75
        elif hit_info['temp_description']:
            annotation_score = 0.5
        else:
            annotation_score = 0.25
        return annotation_score

    def get_intersecting_sources(self,combo,query_hits):
        res=1
        for hit in combo:
            for q_hit in query_hits:
                hit_source,hit_name,hit_info=hit
                q_hit_source,q_hit_name,q_hit_info=q_hit
                if hit_source!= q_hit_source:
                    if self.is_hit_match(hit_info,hit_name,hit_source,q_hit_info,q_hit_name,q_hit_source):
                        res+=1
        return res



    def get_combo_score_Consensus(self,combo,query_length,minX_evalue,maxX_evalue,hmm_sources):
        '''
        metrics for consensus:
        intersecting sources = amount of different sources that intersect with the combo
        combo coverage  = percentage of the sequence covered by the combo
        combo evalue = summed scaled evalue of the combo. Is averaged
        combo weight = summed hmm weight of the combo. Better sources should give better combinations of hits . Is averaged
        combo annotation score = score of the annotation (0.25 if just the hmm name, 0.5 if just the description, 0.75 if just identifiers, 1 if identifiers and description
        '''
        total_evalue=0
        total_coverage=0
        total_weight=0
        annotation_score=0
        for hit in combo:
            hmm_file, hmm_hit, hit_info = hit
            coverage = hit_info['query_end'] - hit_info['query_start']+1
            annotation_score+=self.get_annotation_score(hit_info)
            hmm_weight=self.get_hmm_weight(hmm_file)
            if hit_info['evalue']: evalue = self.min_max_scale_normal(log10(hit_info['evalue']),minX_evalue,maxX_evalue)+0.01
            else:               evalue = 2
            total_coverage+=coverage
            total_evalue+=evalue
            total_weight+=hmm_weight
        sources_score=(hmm_sources + annotation_score + total_weight)/3
        combo_score = sources_score * total_evalue * (total_coverage / query_length)
        combo_score /= len(combo)
        return combo_score


    def get_best_hits_Consensus(self, query_hits,query_length,time_limit=60):
        cython_hits,conversion_dict=self.query_hits_to_cython_Consensus(query_hits)
        cython_possible_combos=get_non_overlapping_hits(cython_hits,time_limit=time_limit)
        #sometimes this calculation is not feasible (when we have too many small hits and the query sequence is too long, the calculation of conbinations would take forever - limit of 5mins)
        if not cython_possible_combos: return None
        best_hit_score=0
        best_combo = None
        minX_evalue,maxX_evalue= self.get_min_max_Consensus(cython_possible_combos,conversion_dict)
        for len_combo in cython_possible_combos:
            if len_combo:
                for cython_combo in cython_possible_combos[len_combo]:
                    combo=self.cython_to_query_hits_Consensus(cython_combo,conversion_dict)
                    #we benefit the combo which intersect with other sources
                    hmm_sources = self.get_intersecting_sources(combo,query_hits)
                    hit_score=self.get_combo_score_Consensus(combo,query_length,minX_evalue,maxX_evalue,hmm_sources)
                    if not best_combo:
                        best_combo = combo
                        best_hit_score = hit_score
                    elif hit_score > best_hit_score :
                        best_hit_score = hit_score
                        best_combo = combo
        return best_combo


    def sort_by_evalue_Consensus(self,query_hits):
        return sorted(query_hits, key=lambda k: k[2]['evalue'])

    def is_overlap_Consensus(self, temp_queries,  current_query):
        if not temp_queries or not current_query: return False
        y = range(current_query[2]['query_start'], current_query[2]['query_end'] + 1)
        for t in temp_queries:
            x = range(t[2]['query_start'], t[2]['query_end'] + 1)
            xs = set(x)
            res = xs.intersection(y)
            if res: return True
            if t[1]==current_query[1]:  return True
        return False

    def get_best_hits_approximation_Consensus(self, query_hits):
        query_hits = self.sort_by_evalue_Consensus(query_hits)
        combo=[]
        while query_hits:
            next_hit= query_hits.pop(0)
            if not self.is_overlap_Consensus(combo,next_hit):
                combo.append(next_hit)
        return combo

    #@timeit_function
    def add_to_best_hits(self,best_hits,query_dict):
        hits_merged=set()
        hmm_files_consensus=set()
        hmm_names_consensus=set()
        for best_hit in best_hits:
            best_hit_file,best_hit_name,best_hit_info=best_hit
            hits_merged.add(best_hit_name)
            hmm_files_consensus.add(best_hit_file)
            hmm_names_consensus.add(best_hit_name)
            temp_best_hit_info=dict(best_hit_info)
            first_check=False
            #iterations might change already iterated over best_hit_info, this way we check if there is any info added, if so, we repeat the cycle one more  otherwise we return
            while temp_best_hit_info!=best_hit_info or not first_check:
                for hmm_file in query_dict:
                    for hit_to_test in query_dict[hmm_file]:
                        if self.is_hit_match(query_dict[hmm_file][hit_to_test],hit_to_test,hmm_file, best_hit_info,best_hit_name,best_hit_file):
                            self.add_to_hit(query_dict[hmm_file][hit_to_test],best_hit_info)
                            hits_merged.add(hit_to_test)
                            hmm_files_consensus.add(hmm_file)
                            hmm_names_consensus.add(hit_to_test)
                first_check=True
                temp_best_hit_info = dict(best_hit_info)
        #checking how many hits we managed to merge out of all the hits - coverage_consensus
        total_hits=0
        for hmm_file in query_dict:
            for _ in query_dict[hmm_file]:
                total_hits+=1
        return str(len(hits_merged)),str(total_hits),hmm_files_consensus,hmm_names_consensus

    def is_nlp_match(self,hit1_info_description,hit2_info_description):
        for hit1_d in hit1_info_description:
            for hit2_d in hit2_info_description:
                score = self.get_similarity_score(hit1_d, hit2_d)
                #print('score',score,'best_hit_d',best_hit_d,'hit_to_test_d',hit_to_test_d)
                if score > self.nlp_threshold:
                    return True
        return False

    def pfam_in_description(self,hit1_identifiers,hit1_description,hit2_identifiers,hit2_description):
        pfam_ids1={i.split(':')[1] for i in hit1_identifiers if 'pfam:' in i}
        pfam_ids2={i.split(':')[1] for i in hit2_identifiers if 'pfam:' in i}
        for hit1_d in hit1_description:
            for hit2_i in pfam_ids2:
                if hit2_i in hit1_d: return True
        for hit2_d in hit2_description:
            for hit1_i in pfam_ids1:
                if hit1_i in hit2_d: return True
        return False

    #eggnog is not specific enough with the go identifiers, so if there is more than one, we dont use it for the consensus
    def remove_unspecific_identifiers(self,identifiers,specificity_threshold=5):
        res = set()
        dict_counter={}
        for i in identifiers:
            link_type=i.split(':')[0]
            if link_type not in dict_counter: dict_counter[link_type]=0
            dict_counter[link_type]+=1
        for i in identifiers:
            link_type=i.split(':')[0]
            if dict_counter[link_type]<=specificity_threshold:
                res.add(i)
        return res

    #this makes mantis more efficient since it saves matches in memory, however it also makes it more ram intensive
    def is_hit_match(self,hit1_info,hit1_name,hit1_source,hit2_info,hit2_name,hit2_source):
        #cleaning up memory
        if float(psutil.virtual_memory().percent)>95:self.query_match_hits={}
        #to avoid adding duplicate entries
        sorted_hits=sorted([hit1_source,hit2_source])
        if hit1_source==sorted_hits[0]:
            first_hit_source=hit1_source
            first_hit_name=hit1_name
            second_hit_source=hit2_source
            second_hit_name=hit2_name
        else:
            first_hit_source=hit2_source
            first_hit_name=hit2_name
            second_hit_source=hit1_source
            second_hit_name=hit1_name
        dict_key=first_hit_source+'_'+first_hit_name+'__'+second_hit_source+'_'+second_hit_name
        if dict_key in self.query_match_hits:
            return self.query_match_hits[dict_key]
        if hit1_source==hit2_source and hit1_name == hit2_name:
            self.query_match_hits[dict_key]=False
            return False
        #now to actually check
        hit1_info_description = hit1_info['temp_description']
        hit1_info_identifiers = hit1_info['temp_identifiers']
        hit2_info_description = hit2_info['temp_description']
        hit2_info_identifiers = hit2_info['temp_identifiers']
        temp_hit1_info_identifiers = self.remove_unspecific_identifiers(hit1_info_identifiers)
        temp_hit2_info_identifiers = self.remove_unspecific_identifiers(hit2_info_identifiers)
        if temp_hit1_info_identifiers.intersection(temp_hit2_info_identifiers):
            self.query_match_hits[dict_key]=True
            return True

        if self.pfam_in_description(temp_hit1_info_identifiers,hit1_info_description,temp_hit2_info_identifiers,hit2_info_description):
            self.query_match_hits[dict_key]=True
            return True

        if self.is_nlp_match(hit1_info_description,hit2_info_description):
            self.query_match_hits[dict_key]=True
            return True
        self.query_match_hits[dict_key] = False
        return False



    def add_to_hit(self,hit_to_test_info,best_hit_info):
        best_hit_info['identifiers'].update(hit_to_test_info['identifiers'])
        best_hit_info['description'].update(hit_to_test_info['description'])

    def generate_consensus_line(self,query,query_dict,is_essential,consensus_hits,total_hits,hmm_files_consensus,hmm_names_consensus):
        hmm_hits=';'.join(hmm_names_consensus)
        hmm_files=';'.join(hmm_files_consensus)
        #consensus_coverage is a str with consensus sources/all sources, better as a str instead of float as its easier to understand
        row_start = [query,hmm_files,hmm_hits,consensus_hits,total_hits, '|']
        row_start=[str(i) for i in row_start]
        all_descriptions=set()
        all_identifiers=set()
        for hmm_file in hmm_files_consensus:
            for hmm_hit_name in hmm_names_consensus:
                if hmm_hit_name in query_dict[hmm_file]:
                    description=query_dict[hmm_file][hmm_hit_name]['description']
                    identifiers=query_dict[hmm_file][hmm_hit_name]['identifiers']
                    all_identifiers.update(identifiers)
                    all_descriptions.update(description)
        res = list(row_start)
        sorted_identifiers = sorted(all_identifiers)
        for link in sorted_identifiers:
            if 'enzyme_ec' in link:
                if link not in res:
                    res.append(link)
        # so that description always comes in the end
        for link in sorted_identifiers:
            if 'kegg_map_lineage' not in link:
                if link not in res:
                    res.append(link)
        for link in sorted_identifiers:
            if 'kegg_map_lineage' in link:
                if link not in res:
                    res.append(link)
        if is_essential:
            res.append('is_essential_gene:True')
        for link in all_descriptions:
            if 'description:'+link not in res:
                res.append('description:'+link)
        return res

    def generate_consensus_output(self,interpreted_annotation_tsv,consensus_annotation_tsv,stdout_file_path=None):
        first_line=['Query',
                    'HMM_Files',
                    'HMM_Hits',
                    'Consensus_hits',
                    'Total_hits',
                    #'Ratio_difference_to_maximum_hit_evalue',
                    #'Maximum_group_evalue',
                    #'Group_Coverage',
                    #'Has_Identifiers',
                    #'Description_significance',
                    '|',
                    'Links',
                    ]
        with open(consensus_annotation_tsv,'w+') as file:
            first_line='\t'.join(first_line)
            file.write(first_line+'\n')
            consensus_annotation = self.generate_consensus(interpreted_annotation_tsv,stdout_file_path=stdout_file_path)
            for consensus_query in consensus_annotation:
                line = '\t'.join(consensus_query)
                file.write(line + '\n')

if __name__ == '__main__':

    f='/home/pedroq/Desktop/test_consensus/exact_int.tsv'
    f2='/home/pedroq/Desktop/test_consensus/consensus_annotation.tsv'
    go_terms_path='/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/Resources/Gene_Ontology/go.obo'
    reference_file='/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/Resources/Uniprot/uniprot-filtered-reviewed_yes.tab'
    m=MANTIS_Consensus()
    m.go_terms_path=go_terms_path
    m.uniprot_reference=reference_file
    m.pickled_go_syns = m.go_terms_path + '.pickle_syns'
    m.pickled_go_terms = m.go_terms_path + '.pickle_terms'
    m.pickled_go_dict = m.go_terms_path + '.pickle_dict'
    #m.parse_go_terms()
    m.nlp_threshold=0.8
    m.domain_algorithm='exhaustive'
    m.mantis_hmm_weights={'else':0.7,'tigrfam_merged':1}
    m.n_grams_range=[1]
    m.words_to_remove = ['mainrole','sub1role','protein','proteins',
                                'enzyme','enzymes','putative','activity',
                                'process','unknown','function','functions',
                                'processes'
                                ]
    m.overlap_value=0.1
    m.false_init()



    #a,d=m.generate_consensus_free_text(query_dict)
    #a1=list(a)
    #a1[0].append('hmm4')
    #m.merge_consensus_ids_nlp(a,a)
    #a=m.min_max_scale(log10(4.04908e-71),-16.7,-102.17)
    #m.generate_consensus(f)
    m.generate_consensus_output(f,f2)
    #print(a)