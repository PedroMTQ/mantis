# Additional information

## Reference data
Mantis uses several default data sources, to stop using these data sources add `NA` to the default hmm path you'd like to stop using, for example: `NOGT_hmm_folder=NA`  
The default hmm data sources are:   
- Pfam
- KOfam
- TIGRfam
- Resfams
- EGGNOG

It also uses a list of [essential genes](https://raw.githubusercontent.com/MadsAlbertsen/multi-metagenome/master/R.data.generation/essential.hmm) to identify HMM profile matches that correspond to previously identified essential genes.
To identify the taxonomic lineage it uses NCBI's taxonmic lineage found [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump).
Respective citations can be found in [References and Acknowledgements](https://github.com/PedroMTQ/mantis/wiki/Home#References-and-acknowledgements)


## What is the e-value threshold?

We have found that setting the e-value to 1e-3 produces the best results. While this values gives a bigger list of hits, the most significant hits are still selected.
This can however be set to another value by adding `-et evalue`.
You can set the evalue to 'dynamic' so that the e-value threshold changes according to the query sequence length, such that:
- if query sequence length <= 150:      e-value threshold = 1e<sup>-10</sup>
- if query sequence length >= 250:    e-value threshold = 1e<sup>-25</sup>
- else:                   e-value threshold = 1e<sup>-query sequence length/10</sup>

Since fasta files are split into chunks, the e-value needs to be adjusted. This is why, for the same hmm hit, you will see different e-values from HMMER's domtblout to Mantis.

## How are hits generated?

Possible combination hits are recursively generated similarly to how all possible combinations would be generated.
Combination of hits are chosen according to 3 different metrics, according to what defines a "good combination":
 - High coverage;  
 - Low e-value;  
 - Low number of hits per combinations - many small sized hits are less significant than a few large sized hits;

This is defined by the following equation:

![Combination of hits score equation](Images/score_equation_small.png)

The left fraction's numerator depicts scaled e-value, the denominator the number of hits.  
The right fraction depicts the coverage (percentage of the query sequence covered by our combination) of the combination of hits.  

**This score ranges from 0 to 1, 1 being the highest.**

### Why log10?  
The e-values within our threshold are extremely low, making comparisons between the raw e-values meaningless. By using the log10 we reduce the range of this distribution, allowing comparison. This also avoids the precision loss that is inherent to the arithmetic capacities of all computer systems.

### Why MinxMaxScale?  
Other scaling equations could be applied, but the MinMaxScale transforms our e-values into proportions, which fits well the combination hit score calculation equation used.

![MinMaxScale equation](Images/min_max_scale_equation.png)

The total amount of possible combinations (discarding the empty combination) is 2<sup>N</sup>-1 where N is the amount of hits our query has.
Calculating and evaluating all possible combinations is not computationally feasible:

- 10 hits: 2<sup>10</sup> - 1 = 1023
- 15 hits: 2<sup>15</sup> - 1 = 32767
- 30 hits: 2<sup>30</sup> - 1 = 1073741823

When dealing with queries with a large number of hits, testing (and even generating) all combinations is not feasible.
The adapted algorithm is as follows:
1. Create an initial list of possible combinations - where the combination root is a query hit
2. For each possible combination test possible hits that can be added
3. Can combination be expanded? Go to 2; otherwise go to 4
4. Return possible combinations

Please see the function `add_to_combo` to check code.

However, there are cases when this improved algorithm may not be feasible either (even with Cython):
- Long query sequence
- Low hit coverage
- High number of hits

If the previous algorithm doesn't finish within 1 minute, a lazy algorithm comes into play:
1. Select lowest e-value hit and add it to the combination.
2. Does hit overlap with our combination? If not, add it to combination
3. Are there any hits left? If yes go to 1.

When iterating through the hits, if a hit has the same or a slightly higher e-value but higher coverage, this hit would be considered better than the hit with the lowest e-value; this flexible threshold is defined by the `acceptable_range` variable.
The default value is 0.05 but you can define another value, including 0 to remove this flexibility.

**Important**: keep in mind that the best hits selection only occurs within the same data source HMM (that's why it's important to merge the several HMMs profiles of one data source into one single HMM).  


To see how this works, check `MANTIS_Processor.py`.


## How are hit consensus generated?

As previously explained, Mantis can get non-overlapping hits from the same HMM source. However, since Mantis uses a lot of different HMM sources, there might not be a consensus between them.
That is why it is crucial to find a way to find some sort of consensus between the different HMM sources.  
The initial step is to identify the best hits for the current query using the same algorithm employed for finding non-overlapping hits. This guarantees that the consensus has the best possible hits amongst all HMM sources.  
In the next step all the remaining hits are cycled through and the consensus is extended (given that the hit is similar to the consensus).   

To find the consensus a two-fold approach was used:  
- Consensus between identifiers
- Consensus between free-text description

The consensus between identifiers is simply calculated by identifying intersections between the different sources.
On the other hand, the consensus between the free-text uses a more advanced approach. It uses a text mining approach to compare the free-text between sources and identify possible matches.

The natural language processing of the annotations free-text entails several steps:
1. Annotation free text pre-processing:
    - Split annotation into documents
    - Remove identifiers
    - Standardize punctuation
    - Remove digits that are not attached to a token
    - Standardize ion patterns
    - Replace Roman numerals with Arabic numerals
    - Divide document into groups of tokens
    - Unite certain tokens (for example: “3” should be merged with “polymerase 3”)
2. [Part-of-speech tagging](#part-of-speech-tagging)  
    - pos_tag with [universal tagging](https://explosion.ai/blog/part-of-speech-pos-tagger-in-python) (contextual)
    - [Wordnet](https://wordnet.princeton.edu/) tagging (independent)
    - Choose best tag (Wordnet takes priority)
    - Removal of unwanted tags (determiners, pronouns, particles, and conjunctions)
3. [Token scoring](#token-scoring)
    - Try to find synonyms (wordnet lexicon) shared between the 2 compared documents
    - Build Term frequency- Inverse Document Frequency vectors (TF-IDF)
4. [Similarity analysis](#similarity-analysis)
    - Calculate cosine distance between the two scaled vectors
    - Calculate Jaccard distance between the two sets of identifiers
    - If similarity score is above the 0.8 consider, it a match
5. [Consensus construction](#consensus-construction)

### Part-of-speech tagging

Part-of-speech tagging (POST) is the method of lexically classifying tokens based on their definition and context. In the context of this application, the point is to eliminate tokens that are not relevant to the similarity analysis.  
After pre-processing, tokens are tagged with a custom tagger [SequentialBackOffTagger](https://kite.com/python/docs/nltk.SequentialBackoffTagger) independent of context. This tagger uses [Wordnet’s lexicon](https://wordnet.princeton.edu/) to identify the most common lexical category of any given token.  
Should a token be present in Wordnet’s lexicon, a list of synonyms and their lexical category is generated, for example:  

`[(token,noun),(synonym1,noun) ,(synonym2,verb),(synonym3,adjective),(synonym4,noun)]`

The token is then assigned the most common tag **noun**.  

To adjust this lexicon to biological data, [gene ontology](http://purl.obolibrary.org/obo/go.obo) tokens are also added.  
Untagged tokens are then contextually classified with a [Perceptron tagger](http://wiki.apertium.org/wiki/Perceptron_tagger). The classification obtained from this tagger is not optimal (as a pre-trained classifier is used), however, in the current context this is barely of consequence, as this tagger is merely used as a backup when no other tag is available. Optimally a new model would be trained, but unfortunately this would require heavy time-investment in building a training dataset.  
The tokens tagged as being determiners, pronouns, particles, or conjunctions are removed.

### Token scoring

In this step, tokens are scored based on the “Term frequency- Inverse Document Frequency” technique. This allows the analysis on which tokens are more relevant to a certain annotation, which in turn allows for the identification of other annotations with the same similarly important tokens.

TF-IDF measures the importance of a token to a document in a corpus. To summarize:
- TF - Tokens that appear more often in a document should be more important. This is a local (document wide) metric.
- IDF - tokens that appear in too many documents should be less important. This is a global (corpus wide) metric.

TF-IDF is calculated with the following equation:

![tf_idf_equation](Images/tf_idf_equation.png)

- NT, times token appears in document
- TT, total amount of tokens in document
- TD, total amount of documents
- DT, total amount of times a certain token appears in a document – frequency table

The corpus used to build this metric were all the 561.911 reviewed proteins from Uniprot (as of 2020/04/14). After pre-processing, each protein annotation is split into tokens, and a token frequency table (DT) is calculated and saved into a file.

The TF-IDF score is then locally scaled (min_max scaling relative to the document) so that we can better understand which tokens are more relevant within the analysed document.


### Similarity analysis

Finally, we can then compare annotations from different sources, by calculating the [cosine distance](https://en.wikipedia.org/wiki/Cosine_similarity) between each pair of TF-IDF scaled vectors. Should the tokens they contain and their importance within the document be around the same, the annotations are classified as “identical”.
Identifiers within the free-text description are also taken into account, via the [Jaccard distance metric](https://en.wikipedia.org/wiki/Jaccard_index). A simple intersection is not used as more general identifiers might lead to too many false positives.

### Consensus construction

In this manner we are able to construct groups of hits (from different sources) that match between either via identifiers or free-text descriptions. We then evaluate the quality of each group of consensuses and select the best one, taking into account:  
- Percentage of the sequence covered by the hits in the consensus
- Significance of the hits (e-value) in the consensus
- Significance of the reference datasets
- Number of different reference datasets in the consensus


To see how this works, check `MANTIS_Consensus.py` and `MANTIS_NLP.py`.  



## Notes on efficiency

Mantis was built to be highly efficient and scalable, while at the same time not relying on heuristic techniques. There are several factors that can cause execution bottlenecks, some of which can be circumvented.  
What are the efficiency bottlenecks?  
- [Workflow design](#workflow-design) – the workflow used is the main contributor to poor performance, an iterative approach is simple to implement but substantially inferior to a parallelized implementation.
- **Number of cores** - the number of cores influences speed performance the most, speed scales linearly with the number of cores, this depends on the hardware available;
- **Number of processes** - one would expect the ideal number of processes to be the same as the number of cores. Up to a certain point increasing the number of processes leads to speed gains, after around 2-3 times the number of cores there’s no more speed gain, increasing past this point can actually lead to severe efficiency decreases due to the overhead generated.
- [Sample size](#splitting-sample-into-chunks) – a larger search space (more query sequences to compare against reference profiles) translates into a higher runtime. An iterative execution can’t possibly handle metagenomic samples in a high throughput scenario, it is essential to properly scale execution so that we can attain results within a reasonable amount of time.
- [Amount of HMM profiles](#splitting-hmms-into-chunks) – again, a larger search space (more reference profiles to match against query sequences) translate into a higher run time.
- **HMMER threads** - thread generation leads to some overhead; I've found that the sweet spot to be around 5 HMMER threads (master thread + 4).

### Workflow design

Mantis parallelizes its inner sub-tasks via Python’s multiprocessing module. It achieves this by having a global task queue that continuously releases tasks to workers. During each main task of the annotation workflow, a certain number of workers are recruited which will execute and consume all the tasks within the queue. When a worker has finished its job, it will execute another task from the queue, until there are no more tasks to execute. If the queue is well balanced, minimal idle time (time spent waiting for workers to empty the task queue) can be achieved.

### Splitting sample into chunks

HMMER performance is determined by several factors, one of which is the length of the query sequence, which can be overcome by splitting samples into smaller chunks, allowing parallelization. It’s important to fine-tune the size of the chunks, splitting into very small chunks ensures saturation, but also generates more overhead (due to process spawning and context switching), splitting into bigger chunks generates less overhead but won't ensure saturation.  
A naive approach but efficient approach to splitting the sample into smaller chunks would be to use a moving window (think of k-mers as an analogy).This means that should a sample with have 3005 query sequences and we are splitting into chunks with 500 sequences, we’d have 6 chunks with 500 sequences and a seventh chunk with 5 sequences; hardly an optimal distribution.  
A better approach would be to evenly distribute sequences across the chunks as well as load balance by sequence length, thus creating more even chunks. This increased complexity is slower and more hardware consuming though.  
An elegant middle ground solution is then to only apply load balancing when dealing with small samples, when uneven chunks would more heavily affect time efficiency. On the other hand, when dealing with big samples (more than 200k query sequences - think metagenomic samples), the uneven chunks (now generated through the moving window approach) are unnoticeable in the sheer volume of chunks that have to annotated. Chunk size is dynamically calculated.


### Splitting HMMs into chunks

Having a large enough reference dataset is essential for blind protein similarity searches, especially when the purpose is the functional annotation of unknown MAGs. While some sources can be a few hundred MBs, others can be spawn to tens of GBs. This discrepancy would lead to substantially different HMMER runtimes. To overcome this issue, HMMs are split into chunks, allowing large datasets to be searched in parallel instead of iteratively. Since both the chunks and the original file are kept (for posteriority), there is duplication of data.

Sample and HMMs splitting achieve quasi-optimal queue balancing.

### Downsides of multiprocessing

Multiprocessing sounds amazing on paper, and, for specific cases, it can indeed be a wonderful tool! It does however have some downfalls:
- Higher RAM consumption due to process spawning (normally nothing to be worried about)
- Minor calculation errors due to intermediate roundings. An anecdotal example:
    - when splitting the fasta into chunks: fasta with 20 sequences, hit had an i-evalue of 5.98e-52
    - when splitting the fasta into chunks: fasta with 24 sequences, hit had an i-evalue of 5.7e-52
    - when the fasta was not split: hit had an i-evalue of 5.8 e-52  
- Higher disk memory consumption due to the generation of intermediate files and  splitting HMMs into chunks.

Since these are minor issues, for optimal efficiency, multiprocessing is still widely used within Mantis.


#### Using multiple physical nodes?

Mantis can't run on multiple physical nodes since multiprocessing requires data to be shared between processes. This can't be achieved with the current implementation, but could potentially be introduced with an OpenMPI implementation.
