from dataclasses import dataclass, field, make_dataclass, fields
import re
import json
from collections.abc import Iterable

TC_PATTERN = re.compile(r'(?<![A-Za-z])\(TC\s\d\.[A-Z\-](\.(\d+|\-)){1,2}\)')
KO_PATTERN = re.compile(r'(?<![A-Za-z])K\d{4,}')
DUF_PATTERN = re.compile(r'(?<![A-Za-z])(DUF|duf)\d{2,}')
PFAM_PATTERN = re.compile(r'(?<![A-Za-z])((U|u)?PF|pf)\d{3,}')
COG_PATTERN = re.compile(r'(?<![A-Za-z])(COG|cog)\d{3,}')
ARCOG_PATTERN = re.compile(r'(AR|ar)(COG|cog)\d{3,}')
TIGRFAM_PATTERN = re.compile(r'(?<![A-Za-z])(TIGR)\d{3,}')
GO_PATTERN = re.compile(r'(?<![A-Za-z])GO\d{3,}')
EC_PATTERN = re.compile(r'\d(\.(-|\d{1,3}|([a-zA-Z]\d{1,3}))){2,3}')


@dataclass
class BaseMetadataDocument():
    enzyme_ec: set = field(default_factory=set)
    tcdb: set = field(default_factory=set)
    tigrfam: set = field(default_factory=set)
    kegg_ko: set = field(default_factory=set)
    pfam: set = field(default_factory=set)
    cog: set = field(default_factory=set)
    arcog: set = field(default_factory=set)
    go: set = field(default_factory=set)
    description: set = field(default_factory=set)

    @classmethod
    def from_string(cls, metadata_string: str):
        document = cls()
        document.add_by_search(metadata_string)
        return document

    @classmethod
    def from_sql(cls, metadata_types: list[str], metadata_list: list[str]):
        metadata_dict = {}
        for metadata_idx, metadata_type in enumerate(metadata_types):
            metadata = json.loads(metadata_list[metadata_idx])
            metadata_dict[metadata_type] = set(metadata)
        document = cls(**metadata_dict)
        return document

    def to_sql(self, metadata_types: list[str]):
        res = []
        for metadata_type in metadata_types:
            metadata_set = getattr(self, metadata_type)
            res.append(json.dumps(list(metadata_set)))
        return res

    def to_file(self):
        res = []
        sorted_metatadata_types = sorted([metadata_type.name for metadata_type in fields(self)])
        sorted_metatadata_types.remove('description')
        sorted_metatadata_types.append('description')
        for metadata_type in sorted_metatadata_types:
            metadata_set = getattr(self, metadata_type)
            if metadata_set:
                for metadata in metadata_set:
                    res.append(f'{metadata_type}:{metadata}')
        return res

    def find_ecs(self, string_to_search, required_level=3) -> set:
        res = set()
        # greedy match of confounders
        search = re.finditer(EC_PATTERN, string_to_search)
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

    def find_tcdb(self, string_to_search) -> set:
        res = set()
        search = re.finditer(TC_PATTERN, string_to_search)
        for i in search:
            res.add(i.group())
        return res

    def find_ko(self, string_to_search) -> set:
        res = set()
        # I could do upper and lower case but since it's only one letter, it's not very safe...
        search = re.finditer(KO_PATTERN, string_to_search)
        for i in search:
            res.add(i.group())
        return res

    def find_pfam(self, string_to_search) -> set:
        res = set()

        for pattern in [DUF_PATTERN, PFAM_PATTERN]:
            search = re.finditer(pattern, string_to_search)
            for i in search:
                res.add(i.group())
        return res

    def find_cog(self, string_to_search) -> set:
        res = set()
        search = re.finditer(COG_PATTERN, string_to_search)
        for i in search:
            res.add(i.group())
        return res

    def find_arcog(self, string_to_search) -> set:
        res = set()
        search = re.finditer(ARCOG_PATTERN, string_to_search)
        for i in search:
            res.add(i.group())
        return res

    def find_tigrfam(self, string_to_search) -> set:
        res = set()
        search = re.finditer(TIGRFAM_PATTERN, string_to_search)
        for i in search:
            res.add(i.group())
        return res

    def find_go(self, string_to_search) -> set:
        res = set()
        search = re.finditer(GO_PATTERN, string_to_search)
        for i in search:
            res.add(i.group())
        return res

    def add(self, metadata_type: str, metadata: str | Iterable):
        metadata_variable = getattr(self, metadata_type)
        if isinstance(metadata, str):
            metadata_variable.add(metadata)
        elif isinstance(metadata, Iterable):
            metadata_variable.update(metadata)
        else:
            raise Exception(f'Invalid metadata type: {metadata} (type: {type(metadata)})')

    def add_by_search(self, metadata_string: str):
        if not metadata_string:
            return
        metadata_function_mapping = {
            'enzyme_ec': self.find_ecs,
            'tcdb': self.find_tcdb,
            'tigrfam': self.find_tigrfam,
            'kegg_ko': self.find_ko,
            'pfam': self.find_pfam,
            'cog': self.find_cog,
            'arcog': self.find_arcog,
            'go': self.find_go,
            }
        for metadata_type, metadata_function in metadata_function_mapping.items():
            metadata = metadata_function(metadata_string)
            metadata_variable = getattr(self, metadata_type)
            metadata_variable.update(metadata)


def get_extended_metadata_document(metadata_types: set[str]):
    extra_fields = [(metadata_type, set, field(default_factory=set)) for metadata_type in metadata_types]
    return make_dataclass(
        cls_name='MetadataDocument',
        fields=extra_fields,
        bases=(BaseMetadataDocument,),
    )


def merge_metadata_documents(metadata_documents: list[BaseMetadataDocument]):
    metadata_types = set()
    for metadata_document in metadata_documents:
        metadata_types.update(metadata_type.name for metadata_type in fields(metadata_document))
    new_metadata_document: BaseMetadataDocument = get_extended_metadata_document(metadata_types=metadata_types)()
    for metadata_document in metadata_documents:
        for metadata_type in fields(metadata_document):
            metadata = getattr(metadata_document, metadata_type.name)
            new_metadata_document.add(metadata_type=metadata_type.name,
                                      metadata=metadata)
    return new_metadata_document


if __name__ == '__main__':
    document = BaseMetadataDocument(enzyme_ec={'1.2.3.4'})
    print(document)
    document = BaseMetadataDocument.from_string('This string contains an enzyme EC 1.2.3.4')
    print(document)
    extended_document_factory: BaseMetadataDocument = get_extended_metadata_document(metadata_types=['type_1', 'type_2'])
    print(extended_document_factory)
    document = extended_document_factory.from_string('This string contains an enzyme EC 1.2.3.4')
    print(document)
    print('json', document.to_sql(['enzyme_ec']))
    print('doc', extended_document_factory.from_sql(metadata_types=['type_1', 'enzyme_ec'],
                                                    metadata_list=['["test"]', '["1.2.3.4"]']))
    document_1 = BaseMetadataDocument(enzyme_ec={'1.2.3.4'})
    extended_document_factory: BaseMetadataDocument = get_extended_metadata_document(metadata_types=['type_1', 'type_2'])
    document_2 = extended_document_factory(type_1='type_1_test')
    merged_document = merge_metadata_documents([document_1, document_2])
    print(document_1)
    print(document_2)
    print(merged_document)
