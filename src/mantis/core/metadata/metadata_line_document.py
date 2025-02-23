from dataclasses import dataclass, field, asdict, fields
from mantis.core.metadata.metadata_document import BaseMetadataDocument, merge_metadata_documents, get_extended_metadata_document
import json


def set_converter(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")

@dataclass
class MetadataLineDocument():
    idx: int
    query: str
    reference_file: str
    reference_id: str
    reference_id_accession: str
    evalue: str
    bitscore: str
    direction: str
    query_length: str
    query_start: str
    query_end: str
    reference_start: str
    reference_end: str
    reference_length: str
    is_essential_gene: bool=field(default_factory=list)
    metadata_documents: list[BaseMetadataDocument]=field(default_factory=list)
    metadata: BaseMetadataDocument=field(default=None)

    def __post_init__(self):
        if self.reference_id_accession == '-':
            self.reference_id_accession = None

    def merge_metadata(self):
        if self.metadata:
            return
        self.metadata = merge_metadata_documents(metadata_documents=self.metadata_documents)
        self.metadata_documents = []

    @classmethod
    def from_json(cls, json_line: str):
        document_dict: dict = json.loads(json_line)
        metadata = document_dict.get('metadata')
        if metadata:
            metadata_document = get_extended_metadata_document(metadata_types=metadata.keys())(**metadata)
            document_dict['metadata'] = None
            document_dict['metadata_documents'] = [metadata_document]
        return cls(**document_dict)

    @staticmethod
    def get_headers():
        headers = [header.name for header in fields(MetadataLineDocument)]
        headers.remove('metadata_documents')
        return headers

    def to_json(self):
        self.merge_metadata()
        document_dict = asdict(self)
        document_dict.pop('metadata_documents')
        return json.dumps(document_dict, default=set_converter)



    def to_file(self):
        self.merge_metadata()
        document_dict = asdict(self)
        document_dict['metadata'] = self.metadata.to_file()
        headers = self.get_headers()
        return {k:v for k,v in document_dict.items() if k in headers}


    # TODO finish this. Do we actually need this?
    # TODO @Oskar can you finalize this? I dont know if we really need it
    def to_gff(self):
        self.merge_metadata()
        # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        # verified with http://genometools.org/cgi-bin/gff3validator.cgi
        line_split = integrated_line.index('|')
        line_data, annotation_data = integrated_line[:line_split], integrated_line[line_split + 1:]
        query, ref_file, hit, hit_accession, evalue, bitscore, direction, query_len, query_start, query_end, ref_start, ref_end, ref_len = line_data
        # attributes of gff line:

        if hit_accession != '-':
            attributes = f'Name={self.reference_id};Target={self.reference_id} {self.reference_start} {self.reference_end};Alias={self.reference_id_accession}'
        else:
            attributes = f'Name={self.reference_id};Target={self.reference_id} {self.reference_start} {self.reference_end}'
        notes = f'Note=ref_file:{self.reference_file},ref_len:{self.reference_length}'
        dbxref = []
        ontology_terms = []
        descriptions = []

        for i in annotation_data:
            if not i.startswith('go:'):
                dbxref.append(i)
            elif i.startswith('description:'):
                descriptions.append(i)
            else:
                ontology_terms.append(i)
        if descriptions:
            notes += ',' + ','.join(descriptions)
        if dbxref and ontology_terms:
            # annotations_str = f'Dbxref={",".join(self.metadata.go)};Ontology_term={",".join(self.)}'
            annotations_str = 'Dbxref=' + ','.join(dbxref) + ';' + 'Ontology_term=' + ','.join(ontology_terms)
        elif dbxref and not ontology_terms:
            all_annotations = 'Dbxref=' + ','.join(dbxref)
        elif not dbxref and ontology_terms:
            all_annotations = 'Ontology_term=' + ','.join(ontology_terms)
        gff_line = '\t'.join([query, 'Mantis', 'CDS', query_start, query_end, evalue, direction, '0']) + '\t'
        if all_annotations:
            gff_line += ';'.join([attributes, notes, all_annotations])
        else:
            gff_line += ';'.join([attributes, notes])
        sequence_region_line = f'##sequence-region {query} 1 {query_len}'
        return query, sequence_region_line, gff_line



if __name__ == '__main__':
    document_1 = BaseMetadataDocument(enzyme_ec={'1.2.5.6'})
    document_2 = BaseMetadataDocument.from_string('This string contains an enzyme EC 1.2.3.4')
    file_line = [0,
                 'query',
               'reference_file',
               'reference_id',
               'reference_id_accession',
               'evalue',
               'bitscore',
               'direction',
               'query_start',
               'query_end',
               'query_length',
               'reference_start',
               'reference_end',
               'reference_length']

    line = MetadataLineDocument(*file_line)
    line.metadata_documents.append(document_1)
    line.metadata_documents.append(document_2)
    line.merge_metadata()
    print('merged metadata', line.metadata)
    print('to_file', line.to_file())
    json_line = line.to_json()
    print(MetadataLineDocument.from_json(json_line))
    print(MetadataLineDocument.get_headers())