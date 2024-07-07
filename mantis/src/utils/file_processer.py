import gzip
import os
import pathlib
import shutil
from zipfile import ZipFile

from mantis.src.utils.logger import logger


def concat_files(output_file, list_file_paths):
    logger.info(f'Concatenating files into {output_file}')
    with open(output_file, 'wb') as wfd:
        for f in list_file_paths:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
            # forcing disk write
            wfd.flush()
            os.fsync(wfd.fileno())


def uncompress_archive(source_filepath, extract_path=None, block_size=65536, remove_source=False):
    file_name = os.path.basename(source_filepath)
    dir_path = os.path.dirname(source_filepath)
    if not extract_path:
        extract_path = dir_path
    if '.tar' in file_name:
        unpack_archive(source_file=source_filepath,
                       extract_dir=extract_path,
                       remove_source=remove_source)
    # only for files
    elif '.gz' in file_name:
        gunzip(source_filepath=source_filepath,
               dest_filepath=extract_path,
               block_size=block_size,
               remove_source=remove_source,)
    elif '.zip' in file_name:
        unzip_archive(source_file=source_filepath,
                      extract_dir=extract_path,
                      remove_source=remove_source)
    else:
        logger.error(f'Cannot uncompress file, due to incorrect file extension {source_filepath}')


# this unzips to the same directory!
def gunzip(source_filepath, dest_filepath=None, block_size=65536, remove_source=False):
    if not dest_filepath:
        dest_filepath = source_filepath.strip('.gz')
    if os.path.isdir(dest_filepath):
        file_name = os.path.basename(source_filepath)
        file_name = pathlib.Path(file_name).stem
        dest_filepath = os.path.join(dest_filepath, file_name)
    logger.info(f'Gunzipping {source_filepath} to {dest_filepath}')
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)
        d_file.write(block)
    if remove_source: os.remove(source_filepath)


def unpack_archive(source_file, extract_dir, remove_source=False):
    logger.info(f'Unpacking {source_file} to {extract_dir}')
    shutil.unpack_archive(source_file, extract_dir=extract_dir)
    if remove_source:
        os.remove(source_file)


def unzip_archive(source_file, extract_dir, remove_source=False):
    logger.info(f'Unzipping {source_file} to {extract_dir}')
    with ZipFile(source_file, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)
    if remove_source:
        os.remove(source_file)



def move_file(source_file, dest_file):
    if not os.path.isdir(dest_file):
        if os.path.exists(dest_file):
            os.remove(dest_file)
    try:
        os.rename(source_file, dest_file)
    except Exception:
        shutil.move(source_file, dest_file)


def copy_file(source_file, dest_file):
    if not os.path.isdir(dest_file):
        if os.path.exists(dest_file):
            os.remove(dest_file)
    shutil.copyfile(source_file, dest_file)


def remove_file(source_file):
    if os.path.exists(source_file):
        os.remove(source_file)
