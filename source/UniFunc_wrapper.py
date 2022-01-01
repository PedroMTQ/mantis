try:
    from source.utils import unifunc_downloaded,download_unifunc
except:
    from utils import unifunc_downloaded,download_unifunc

try:
    if not unifunc_downloaded():
        download_unifunc()
    from Resources.UniFunc.UniFunc import UniFunc
except:
    from UniFunc.UniFunc import UniFunc


def test_nlp():
    nlp = UniFunc_wrapper()
    str1 = 'Responsible for trypanothione reduction'
    str2 = 'Protein associated with trypanothione reductase activity'
    nlp.get_similarity_score(str1, str2, verbose=True)

#this is basically a wrapper for UniFunc
class UniFunc_wrapper(UniFunc):
    def __init__(self):
        UniFunc.__init__(self)

if __name__ == '__main__':
    nlp = UniFunc_wrapper()
