from pathlib import Path
import csv
from mantis.core.metadata.line_document import MetadataLineDocument


class MetadataFileWriter():
    JSON_EXTENSIONS = ['.jsonl', '.json']
    CSV_EXTENSIONS = ['.csv', '.tsv']
    def __init__(self, file_path):
        self.file = open(file_path, 'w+')
        self.file_extension = Path(file_path).suffix
        self.__writer = None
        if self.file_extension in self.CSV_EXTENSIONS:
            if self.file_extension == '.tsv':
                delimiter = '\t'
            elif self.file_extension == '.csv':
                delimiter = ','
            self.__writer = csv.DictWriter(self.file,
                                           fieldnames=MetadataLineDocument.get_headers(),
                                           delimiter=delimiter)
            self.__writer.writeheader()
        elif self.file_extension in self.JSON_EXTENSIONS:
            pass
        else:
            raise Exception(f'Invalid file extension: {self.file_extension}')

    def write(self, line_document: MetadataLineDocument):
        if self.file_extension in self.JSON_EXTENSIONS:
            self.file.write(f'{line_document.to_json()}\n')
        elif self.file_extension in self.CSV_EXTENSIONS:
            self.__writer.writerow(line_document.to_file())

    def __del__(self):
        self.file.close()


class MetadataGffFileWriter():
    def __init__(self, file_path):
        gff_file_path = f'{Path(file_path).stem}.gff'
        self.file = open(gff_file_path, 'w+')
        self.file.write('##gff-version 3' + '\n')


    def write(self, line_document: MetadataLineDocument):
        self.file.write(f'{line_document.to_gff()}\n')


    def __del__(self):
        self.file.close()
