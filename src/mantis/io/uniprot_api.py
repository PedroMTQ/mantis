
import re
from time import sleep
from urllib.parse import parse_qs, urlencode, urlparse

import requests
from requests.adapters import HTTPAdapter, Retry


class UnitprotApi():
    '''
    This module was copied from https://www.uniprot.org/help/id_mapping and adapted to our needs
    '''
    API_URL = "https://rest.uniprot.org"
    POLLING_INTERVAL = 3
    def __init__(self) -> None:
        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=retries))

    def get_next_link(self, headers: list):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)


    def check_id_mapping_results_ready(self, job_id):
        c = 0
        url = f"{self.API_URL}/idmapping/status/{job_id}"
        while True:
            r = self.session.get(url)
            r.raise_for_status()
            j = r.json()
            if "jobStatus" in j:
                if j["jobStatus"] == "RUNNING":
                    print(f"Retrying in {self.POLLING_INTERVAL}s")
                    sleep(self.POLLING_INTERVAL)
                else:
                    raise Exception(r["jobStatus"])
            else:
                return bool(j["results"] or j["failedIds"])
            if c == 10:
                raise ConnectionError(url)
            c += 1


    def get_batch(self, batch_response, file_format):
        batch_url = self.get_next_link(batch_response.headers)
        while batch_url:
            batch_response = self.session.get(batch_url)
            batch_response.raise_for_status()
            yield self.decode_results(batch_response, file_format)
            batch_url = self.get_next_link(batch_response.headers)


    def combine_batches(self, all_results, batch_results, file_format):
        if file_format == "json":
            for key in ("results", "failedIds"):
                if key in batch_results:
                    if batch_results[key]:
                        all_results[key] += batch_results[key]
        else:
            return all_results + batch_results
        return all_results


    def get_id_mapping_results_link(self, job_id):
        url = f"{self.API_URL}/idmapping/details/{job_id}"
        r = self.session.get(url)
        r.raise_for_status()
        return r.json()["redirectURL"]

    def decode_results(self, response, file_format):
        if file_format == "json":
            return response.json()
        elif file_format == "tsv":
            return [line for line in response.text.split("\n") if line]
        return response.text

    def get_id_mapping_results_search(self, url: str):
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
        r = self.session.get(url)
        r.raise_for_status()
        results = self.decode_results(r, file_format)
        for batch in self.get_batch(r, file_format):
            results = self.combine_batches(results, batch, file_format)
        return results

    def submit_id_mapping(self, from_db, to_db, ids):
        r = requests.post(url=f"{self.API_URL}/idmapping/run",
                          data={"from": from_db,
                                 "to": to_db,
                                   "ids": ",".join(ids)},
                                   )
        r.raise_for_status()
        return r.json()["jobId"]





