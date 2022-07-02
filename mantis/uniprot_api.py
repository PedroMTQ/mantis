'''
This module was copied from https://www.uniprot.org/help/id_mapping
'''

import re
from time import sleep
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

try:
    from mantis.exceptions import *
    from mantis.utils import *
except:
    from exceptions import *
    from utils import *

POLLING_INTERVAL = 3

API_URL = "https://rest.uniprot.org"

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def submit_id_mapping(from_db, to_db, ids):
    r = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    r.raise_for_status()
    return r.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    c = 0
    while True:
        r = session.get(f"{API_URL}/idmapping/status/{job_id}")
        r.raise_for_status()
        j = r.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                sleep(POLLING_INTERVAL)
            else:
                raise Exception(r["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])
        if c == 10:
            kill_switch(ConnectionError, url)
        c += 1


def get_batch(batch_response, file_format):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results:
                if batch_results[key]:
                    all_results[key] += batch_results[key]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    r = session.get(url)
    r.raise_for_status()
    return r.json()["redirectURL"]


def decode_results(response, file_format):
    if file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    return response.text


def print_progress_batches(batch_index, size, total):
    n = min((batch_index + 1) * size, total)
    print(f"Fetched: {n} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    r = session.get(url)
    r.raise_for_status()
    results = decode_results(r, file_format)
    for i, batch in enumerate(get_batch(r, file_format)):
        results = combine_batches(results, batch, file_format)
    return results
