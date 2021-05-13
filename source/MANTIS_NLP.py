try:
    from source.utils import unifunc_downloaded,download_unifunc
except:
    from utils import unifunc_downloaded,download_unifunc

try:
    if not unifunc_downloaded():
        download_unifunc()
    from Resources.UniFunc.source import UniFunc
except:
    from UniFunc.source import UniFunc


def test_nlp():
    nlp = MANTIS_NLP()
    str1 = 'Responsible for trypanothione reduction'
    str2 = 'Protein associated with trypanothione reductase activity'
    nlp.get_similarity_score(str1, str2, verbose=True)

#this is basically a wrapper for UniFunc
class MANTIS_NLP(UniFunc):
    def __init__(self):
        UniFunc.__init__(self)

if __name__ == '__main__':
    nlp = MANTIS_NLP()
