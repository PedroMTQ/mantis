try:
    from mantis.utils import unifunc_downloaded, download_unifunc
except:
    from utils import unifunc_downloaded, download_unifunc

try:
    if not unifunc_downloaded():
        download_unifunc()
    from Resources.UniFunc.unifunc.source import UniFunc
except:
    from UniFunc.unifunc.source import UniFunc


def test_nlp():
    nlp = UniFunc_wrapper()
    str1 = 'Responsible for trypanothione reduction'
    str2 = 'Protein associated with trypanothione reductase activity'
    res = nlp.get_similarity_score(str1, str2, verbose=True)
    print(res)


# this is basically a wrapper for UniFunc
class UniFunc_wrapper(UniFunc):
    def __init__(self):
        UniFunc.__init__(self)


if __name__ == '__main__':
    nlp = UniFunc_wrapper()
    test_nlp()
