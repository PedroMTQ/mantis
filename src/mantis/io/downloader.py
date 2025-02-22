import os
import shutil
import urllib.request as request
from contextlib import closing

import requests

from mantis.utils.exceptions import DownloadFailed
from mantis.io.logger import logger


def download_file_http(url, file_path, c, ctx):
    if c > 5:
        download_file_http_failsafe(url, file_path, ctx)
    else:
        if ctx:
            with requests.get(url, stream=True, verify=False) as r:
                with open(file_path, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)
        else:
            with requests.get(url, stream=True) as r:
                with open(file_path, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)


# slower but safer
def download_file_http_failsafe(url, file_path, ctx):
    with requests.Session() as session:
        if ctx:
            session.verify = False
        get = session.get(url, stream=True)
        if get.status_code == 200:
            with open(file_path, 'wb') as f:
                for chunk in get.iter_content(chunk_size=1024):
                    f.write(chunk)


def download_file_ftp(url, file_path, ctx):
    with closing(request.urlopen(url, context=ctx)) as r:
        with open(file_path, 'wb') as f:
            shutil.copyfileobj(r, f)


def download_file(url, output_folder='', retry_limit=10):
    file_path = output_folder + url.split('/')[-1]
    ctx = None
    try:
        target_file = request.urlopen(url)
    except Exception:
        try:
            import ssl
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE
            target_file = request.urlopen(url, context=ctx)
        except Exception:
            logger.error('Cannot download target url', url)
            return
    target_size = target_file.info()['Content-Length']
    transfer_encoding = target_file.info()['Transfer-Encoding']
    if target_size:
        target_size = int(target_size)
    if os.path.exists(file_path):
        if transfer_encoding == 'chunked':
            return
        elif os.stat(file_path).st_size == target_size:
            logger.info(f'Not downloading from {url} since file was already found!')
            return
        else:
            os.remove(file_path)
    logger.info(f'Downloading from {url} . The file will be kept in {output_folder}')
    c = 0
    while c <= retry_limit:
        if 'ftp' in url:
            try:
                download_file_ftp(url, file_path, ctx)
            except Exception:
                try:
                    download_file_http(url, file_path, c, ctx)
                except Exception:
                    pass
        else:
            try:
                download_file_http(url, file_path, c, ctx)
            except Exception:
                pass
        if transfer_encoding == 'chunked':
            return
        if os.path.exists(file_path):
            if os.stat(file_path).st_size == target_size:
                return
        c += 1
    raise DownloadFailed(url)
