from dataclasses import dataclass, field
from mantis.core.metadata.metadata_document import BaseMetadataDocument, merge_metadata_documents



@dataclass
class LineDocument():
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
    is_essential_gene: bool=field(default=False)
    metadata_documents: list[BaseMetadataDocument]=field(default_factory=list)
    metadata: BaseMetadataDocument=field(default=None)

    def __post_init__(self):
        if self.reference_id_accession == '-':
            self.reference_id_accession = None

    def merge_metadata(self):
        self.metadata = merge_metadata_documents(metadata_documents=self.metadata_documents)

    def to_file(self):
        res = [self.query,
               self.reference_file,
               self.reference_id,
               self.reference_id_accession,
               self.evalue,
               self.bitscore,
               self.direction,
               self.query_start,
               self.query_end,
               self.query_length,
               self.reference_start,
               self.reference_end,
               self.reference_length,
               '|']
        res.extend(self.metadata.to_file())
        return res

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

    line = LineDocument(*file_line)
    line.metadata_documents.append(document_1)
    line.metadata_documents.append(document_2)
    line.merge_metadata()
    print('merged metadata', line.metadata)
    print(line.to_file())