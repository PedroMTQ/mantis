import os
import re
import shutil
import subprocess
from datetime import datetime
from functools import wraps
from math import ceil
from pickle import dump as pickle_dump
from pickle import load as pickle_load
from time import sleep, time

import psutil

from mantis.src.settings import (
    AVAILABLE_RAM,
    CYTHON_FOLDER,
    ENVIRONMENT_CORES,
    PROCESS_OVERHEAD,
    PSUTIL_EXCEPTIONS,
    WORKER_PER_CORE,
)
from mantis.settings import BATCH_SIZE


def batch_yielder(data_yielder):
    batch = []
    for i in data_yielder:
        if len(batch) < BATCH_SIZE:
            batch.append(i)
        else:
            yield batch
            batch = []
            batch.append(i)
    if batch:
        yield batch


def kill_switch(error_type, message='', flush=False, file=None):
    import signal
    master_pid = os.getpid()
    if file and message:
        print(message, flush=flush, file=file)
    elif message:
        print(error_type + message)
    else:
        print(error_type)
    sleep(5)
    os.kill(master_pid, signal.SIGKILL)




def estimate_number_workers_setup_database(n_hmms_to_setup,
                                           percentage_allowed_overhead=0.8,
                                           minimum_jobs_per_worker=3,
                                           user_cores=None):
    if user_cores:        return user_cores * WORKER_PER_CORE
    # 0.1GB RAM per hmm setup
    ram_per_setup = 0.1
    maximum_overhead_workers = ceil((AVAILABLE_RAM * percentage_allowed_overhead) / PROCESS_OVERHEAD)
    maximum_workers = ceil(AVAILABLE_RAM / (ram_per_setup + PROCESS_OVERHEAD))
    environment_workers = ENVIRONMENT_CORES * WORKER_PER_CORE
    return ceil(min([environment_workers, n_hmms_to_setup, maximum_workers, maximum_overhead_workers]))


def estimate_number_workers_process_output(n_chunks, searchout_per_chunks=1, user_cores=None):
    if user_cores: return user_cores
    environment_workers = ENVIRONMENT_CORES * WORKER_PER_CORE
    return min([environment_workers, n_chunks * searchout_per_chunks])



def get_seqs_count(target_sample):
    total_seqs = 0
    with open(target_sample) as file:
        line = file.readline()
        while line:
            if line[0] == '>': total_seqs += 1
            line = file.readline()
    return total_seqs

def read_profile(hmm_file):
    res = []
    line = hmm_file.readline()
    while line:
        # print(repr(line),line[0:2])
        res.append(line)
        if '//' in line[0:2]:
            return res
        line = hmm_file.readline()


def get_ref_in_folder(folder_path):
    temp_path = add_slash(folder_path)
    if os.path.exists(temp_path):
        if temp_path[0:2] != 'NA':
            files = os.listdir(temp_path)
            for f in files:
                if f.endswith('.hmm'):
                    return temp_path + f
                if f.endswith('.dmnd'):
                    return temp_path + f
    return None


def get_hmm_profile_count(hmm_path, stdout=None):
    res = 0
    if hmm_path.endswith('.hmm'):
        with open(hmm_path) as hmm_file:
            line = hmm_file.readline()
            while line:
                if 'NAME' in line:
                    res += 1
                    if stdout:
                        if not res % 10000:
                            print('Counted ' + str(res) + ' profiles on ' + hmm_path + ' so far.', flush=True,
                                  file=stdout)
                line = hmm_file.readline()
    return res


def get_chunks_path(hmm_path: str) -> list:
    # this will check if the hmm have been split into chunks, if so it will retrieve the chunks instead
    res = []
    hmm_folder = os.path.dirname(hmm_path)
    if 'chunks' in os.listdir(hmm_folder):
        for chunk_file in os.listdir(os.path.join(hmm_folder, 'chunks')):
            if chunk_file.endswith('.hmm'):
                res.append(os.path.join(hmm_folder, 'chunks', chunk_file))
        return res
    else:
        return [hmm_path]


def round_to_digit(number):
    number_str = str(int(number))
    res = number_str[0:1]
    to_add = 1
    if len(number_str) > 1:
        if int(number_str[1]) >= 5:
            to_add += 1
            res += '5'
    for _ in range(len(number_str) - to_add):
        res += '0'
    return int(res)


def chunk_generator(to_chunk, chunk_size):
    for i in range(0, len(to_chunk), chunk_size):
        yield to_chunk[i:i + chunk_size]


def chunk_generator_load_balanced(list_ordered, chunk_size, time_limit=60):
    '''
    :param list_chunks: list of keys which will correpond to sequences IDs . list should have been ordered by sequence length
    :param time_limit: when dealing with a lot of sequences, this method may take a while, if it takes too long we raise TimeoutError
    #1M sequences take up about 80-90 seconds to create the chunks, which is still kinda okay. This should be below most metagenomic samples
    '''
    n_chunks = ceil(len(list_ordered) / chunk_size)
    res = []
    direction_chunks = {}
    for i in range(n_chunks):
        res.append([])
        direction_chunks[i] = True
    chunk_index = 0
    start = time()
    while list_ordered:
        if time_limit:
            if time() - start > time_limit: raise TimeoutError
        if direction_chunks[chunk_index]:
            chunk_val = list_ordered.pop(0)
            direction_chunks[chunk_index] = False
        else:
            chunk_val = list_ordered.pop(-1)
            direction_chunks[chunk_index] = True
        res[chunk_index].append(chunk_val)
        if chunk_index == n_chunks - 1:
            chunk_index = 0
        else:
            chunk_index += 1
    return res


########Processing protein fasta
def remove_temp_fasta(temp_fasta_path, db):
    if f'missing_annotations.{db}.tmp' in temp_fasta_path:
        os.remove(temp_fasta_path)


def is_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        line = f.readline()
        if line[0] == '>':
            return True
    return False


def get_seq_names(protein_fasta_path):
    res = set()
    with open(protein_fasta_path, 'r') as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                query = line.replace('>', '').strip()
                print(query)
                res.add(query)
            line = file.readline()
    return res


def process_protein_fasta_line(res, query, fasta_line, start_recording):
    # for the first > line
    if '>' in fasta_line and not start_recording:
        fasta_line = fasta_line.replace('\'', '')
        fasta_line = fasta_line.replace('>', '')
        fasta_line = fasta_line.replace('\"', '')
        new_query = fasta_line.split()[0]
        start_recording = True
        res[new_query] = ''
        return res, new_query, start_recording
    # for posterior > lines
    elif '>' in fasta_line and start_recording:
        fasta_line = fasta_line.replace('\'', '')
        fasta_line = fasta_line.replace('>', '')
        fasta_line = fasta_line.replace('\"', '')
        start_recording = True
        new_query = fasta_line.split()[0]
        res[new_query] = ''
        return res, new_query, start_recording
    # to get the sequences
    elif start_recording:
        fasta_line = fasta_line.replace('\"', '')
        fasta_line = fasta_line.replace('\'', '')
        res[query] += fasta_line.strip().upper()
        return res, query, start_recording


# this is only used for internal temporary files
def read_protein_fasta(protein_fasta_path):
    res = {}
    with open(protein_fasta_path, 'r') as file:
        line = file.readline()
        if line[0] != '>': kill_switch(InvalidFastaFormat, 'Please make sure there are no blank lines')
        start_recording = False
        query = None
        while line:
            line = line.strip('\n')
            if line:
                res, query, start_recording = process_protein_fasta_line(res=res,
                                                                         query=query,
                                                                         fasta_line=line,
                                                                         start_recording=start_recording)
            line = file.readline()
    return res


def yield_target_seqs(protein_fasta_path, target_seqs):
    query = None
    seq = []
    with open(protein_fasta_path, 'r') as file:
        line = file.readline()
        if line[0] != '>': kill_switch(InvalidFastaFormat, 'Please make sure there are no blank lines')
        while line:
            if line.startswith('>'):
                if query:
                    if query in target_seqs:
                        seq = ''.join(seq).upper()
                        yield f'>{query}\n{seq}\n'
                    seq = []
                query = line.replace('>', '').strip()
            else:
                seq.append(line.strip())
            line = file.readline()
        if query:
            if query in target_seqs:
                seq = ''.join(seq).upper()
                yield f'>{query}\n{seq}\n'


def yield_target_metadata(metadata_file, target_seqs):
    with open(metadata_file, 'r') as file:
        for line in file:
            seq = line.split('\t')[0]
            if seq in target_seqs:
                yield line


# low memory footprint_version
def read_protein_fasta_generator(protein_fasta_path):
    query = None
    seq = []
    with open(protein_fasta_path, 'r') as file:
        line = file.readline()
        if line[0] != '>': kill_switch(InvalidFastaFormat, 'Please make sure there are no blank lines')
        while line:
            if line.startswith('>'):
                if query:
                    if seq:
                        yield query, ''.join(seq).upper()
                    seq = []
                query = line.replace('>', '').strip()
            else:
                seq.append(line.strip())
            line = file.readline()
        if query:
            if seq:
                yield query, ''.join(seq).upper()


# low memory footprint_version
def write_fasta_generator(seqs_generator, fasta_file):
    with open(fasta_file, 'w+') as file:
        for query, seq in seqs_generator:
            chunks = [seq[x:x + 60] for x in range(0, len(seq), 60)]
            chunk_str = '\n'.join(i for i in chunks)
            file.write('>' + query + '\n' + chunk_str + '\n')


def add_slash(path_folder):
    if not path_folder: return path_folder
    if path_folder[-1] != SPLITTER: return path_folder + SPLITTER
    return path_folder


def get_path_level(path, level=1, remove_extension=False):
    temp = path.strip(SPLITTER)
    res = temp.split(SPLITTER)
    res = res[-level]
    if level == 1 and remove_extension: res = res.split('.')[0]
    return res


def print_cyan(*args, flush=False, file=None):
    for i in args:  print('# ', i, flush=flush, file=file)


def red(text, flush=False, file=None):
    print('\033[31m', text, '\033[0m', flush=flush, sep='', file=file)


def green(text, flush=False, file=None):
    print('\033[32m', text, '\033[0m', flush=flush, sep='', file=file)


def yellow(text, flush=False, file=None):
    print('\033[33m', text, '\033[0m', flush=flush, sep='', file=file)


def timeit_function(f):
    @wraps(f)
    def wrapper(self, *args, **kwargs):
        start_time = time()
        res = f(self, *args, **kwargs)
        print('This function', f.__name__, 'took', time() - start_time, 'seconds to run')
        return res

    return wrapper


def timeit_class(f):
    @wraps(f)
    def wrapper(self, *args, **kwargs):
        res = f(self, *args, **kwargs)
        datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print_cyan('This ' + type(self).__name__ + ' process finished running at ' + datetime_str + ' and took ' + str(
            int(time() - self.start_time)) + ' seconds to complete.', flush=True, file=self.redirect_verbose)
        return res

    return wrapper






def merge_profiles(folder_path, output_file, stdout_file=None):
    print('Merging profiles in ', folder_path, flush=True, file=stdout_file)
    main_folder = get_folder(output_file)
    old_files = os.listdir(main_folder)
    for file in old_files:
        if '.hmm' in file:
            os.remove(main_folder + file)
    list_dir = os.listdir(folder_path)
    profiles = [folder_path + SPLITTER + i for i in list_dir if '.hmm' in i.lower()]
    concat_files(output_file, profiles, stdout_file=stdout_file)
    if os.path.exists(folder_path):      shutil.rmtree(folder_path)


def merge_redundant_profiles(output_file, list_file_paths, stdout_file=None):
    already_added = set()
    current_hmm = ''
    current_hmm_name = ''
    print('Merging profiles to ', output_file, flush=True, file=stdout_file)
    with open(output_file, 'w+') as outfile:
        for f in list_file_paths:
            with open(f) as infile:
                line = infile.readline()
                while line:
                    if line[0:4] == 'NAME':
                        if current_hmm and current_hmm_name not in already_added:
                            already_added.add(current_hmm_name)
                            outfile.write(current_hmm)
                        current_hmm = ''
                        current_hmm_name = line.strip('\n')
                        current_hmm_name = current_hmm_name.replace('NAME', '')
                        current_hmm_name = current_hmm_name.strip()
                    current_hmm += line
                    line = infile.readline()


def merge_redundant_sql_annotations(output_file, list_file_paths, stdout_file=None):
    already_added = set()
    print('Merging sql_annotations to ', output_file, flush=True, file=stdout_file)
    with open(output_file, 'w+') as outfile:
        for f in list_file_paths:
            with open(f) as infile:
                line = infile.readline()
                while line:
                    hmm_name = line.split('\t')[0]
                    if hmm_name not in already_added:
                        already_added.add(hmm_name)
                        outfile.write(line)
                    line = infile.readline()



def get_available_ram_percentage(worker_status, user_memory=None):
    workers_ram = 0
    worker_ids = list(worker_status.keys())
    for worker_id in worker_ids:
        try:
            workers_ram += psutil.Process(worker_id).memory_info().rss
        except PSUTIL_EXCEPTIONS:
            worker_status.pop(worker_id)

    workers_ram /= 1024 ** 3
    if user_memory:
        workers_ram_percent = 100 * workers_ram / user_memory
    else:
        workers_ram_percent = 100 * workers_ram / AVAILABLE_RAM
    return 100 - workers_ram_percent


def get_child_workers(master_pid, wanted_name=None):
    '''
    master_pid is the main mantis process (lvl0) which spawns all other multiprocessing processes (lvl1)
    these processes (lvl1) will then run a X command which in turn are themselves processes (lvl2)
    so what we want are not the lvl1 processes but the lvl2 processes. We can know these by their process_name / wanted_name
    '''
    res = set()
    # if not psutil.pid_exists(master_pid): return res
    try:
        master_process = psutil.Process(pid=master_pid)
        children = master_process.children(recursive=True)
        for process in children:
            try:
                if wanted_name:
                    if wanted_name == process.name():
                        res.add(process.pid)
                else:
                    res.add(process.pid)
            except PSUTIL_EXCEPTIONS:
                pass
    except PSUTIL_EXCEPTIONS:
        return res
    return sorted(res)


def count_running_workers(worker_status):
    res = 0
    for w in worker_status:
        if worker_status[w]: res += 1
    return res


def get_workers_status(workers_ids):
    res = {}
    for worker_pid in workers_ids:
        try:
            process = psutil.Process(pid=worker_pid)
            if process.status() == 'running':
                res[worker_pid] = True
            else:
                res[worker_pid] = False
        except PSUTIL_EXCEPTIONS:
            pass
    return res


def kill_workers(master_pid, worker, worker_pid):
    sleep(5)
    try:
        if psutil.pid_exists(master_pid):
            psProcess = psutil.Process(pid=master_pid)
            print(
                '###### Ran out of memory, please increase available memory. Quitting now! ######',
                flush=True)
            psProcess.kill()
            os.waitpid(master_pid, os.WEXITED)
    except PSUTIL_EXCEPTIONS:
        pass
    try:
        worker.kill()
        os.waitpid(worker_pid, os.WEXITED)
    except PSUTIL_EXCEPTIONS:
        pass


def resume_workers(child_status, n_workers):
    ordered_ids = sorted(child_status.keys())
    if count_running_workers(child_status) == n_workers: return
    for worker_pid in ordered_ids:
        try:
            process = psutil.Process(pid=worker_pid)
            if n_workers > 0:
                process.resume()
                child_status[worker_pid] = True
                n_workers -= 1
            else:
                process.suspend()
                child_status[worker_pid] = False
        except PSUTIL_EXCEPTIONS:
            child_status.pop(worker_pid)


def suspend_workers(child_status, n_workers):
    ordered_ids = sorted(child_status.keys(), reverse=True)
    if not n_workers: return
    if len(child_status) - count_running_workers(child_status) == n_workers: return
    for worker_pid in ordered_ids:
        try:
            process = psutil.Process(pid=worker_pid)
            if n_workers > 0:
                process.suspend()
                child_status[worker_pid] = False
                n_workers -= 1
            else:
                process.resume()
                child_status[worker_pid] = True
        except PSUTIL_EXCEPTIONS:
            child_status.pop(worker_pid)


# in order to not run out of memory when running HMMER we control how the processes are running, suspending them if needed
def run_command_managed(command, master_pid, get_output=False, stdout_file=None, wanted_child=None, user_memory=None,
                        shell=False, join_command=False):
    if join_command:        command = ' '.join(command)
    if get_output:
        # popen launches and doesnt wait for finish
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell)
    elif stdout_file:
        process = subprocess.Popen(command, stdout=stdout_file, stderr=stdout_file, shell=shell)
    else:
        process = subprocess.Popen(command, shell=shell)
    process_pid = process.pid
    psProcess = psutil.Process(pid=process_pid)
    # hmmsearch uses a lot of ram when starting up. we need to avoid suspending them when that happens, otherwise we spike memory due to a lot of suspended workers
    # RAM usage usually starts very low, then a huge spike (e.g. 2-5GB for NOGG) and then down again, we can apply a moving average
    moving_average = 0
    c = 0
    spike = False
    spiked = False
    while psutil.pid_exists(process_pid):
        # MB
        process_memory = psProcess.memory_info().rss / 2 ** 20
        if wanted_child == 'hmmsearch':
            # here is when the first spike occurs
            if moving_average:
                # minimum of 100mb for spike detection
                # we only care about the initial spike
                if moving_average > 100 and not spiked:
                    if (process_memory - moving_average) / moving_average > 10:
                        spike = True
                        spiked = True
                # here when the spike has passed/is about to pass
                if process_memory - moving_average < 0:
                    spike = False
            moving_average = moving_average * c
            c += 1
            moving_average = (moving_average + process_memory) / c
            # print(process_memory,moving_average,spike,spiked,command,flush=True)
        if not psutil.pid_exists(master_pid):
            kill_workers(master_pid, process, process_pid)
        if psProcess.status() == 'zombie':
            process.communicate()
            return process
        child_workers = get_child_workers(master_pid=master_pid, wanted_name=wanted_child)
        child_status = get_workers_status(child_workers)
        available_ram_percentage = get_available_ram_percentage(child_status, user_memory)
        # if we have more than 20% ram available we let it run
        if available_ram_percentage >= 20:
            resume_workers(child_status, n_workers=len(child_status))
        # if we have only 10-20% ram available, we suspend half the processes
        elif available_ram_percentage >= 10 and available_ram_percentage < 20:
            workers_to_resume = ceil(len(child_status) / 2)
            if len(child_status) % 2:
                workers_to_suspend = workers_to_resume - 1
            else:
                workers_to_suspend = workers_to_resume
            resume_workers(child_status, n_workers=workers_to_resume)
            if count_running_workers(child_status) > 1:
                if not spike:
                    suspend_workers(child_status, n_workers=workers_to_suspend)
        # if there's not enough memory to run even a single process, we crash execution
        elif available_ram_percentage < 1:
            kill_workers(master_pid, process, process_pid)
        # if we have less than 10% we suspend all but 1
        else:
            # if we have only 1 running worker we let it run
            if count_running_workers(child_status) < 1:
                resume_workers(child_status, n_workers=1)
            # if we have more than one we suspend all but one
            elif count_running_workers(child_status) > 1:
                if not spike:
                    workers_to_suspend = len(child_status) - 1
                    suspend_workers(child_status, workers_to_suspend)
                    resume_workers(child_status, n_workers=1)
        sleep(1)
    return process


def run_command_simple(command, get_output=False, shell=False, join_command=False):
    if join_command:
        command = ' '.join(command)
    # TODO add console output redirection to output
    if get_output:
        # run launches popen and waits for finish
        process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell)
    else:
        process = subprocess.run(command, shell=shell)
    return process


def run_command(command, get_output=False, master_pid=None, wanted_child=None, user_memory=None,
                shell=False, join_command=False):
    command_list = command.split()
    if master_pid:
        return run_command_managed(command=command_list, get_output=get_output,
                                   master_pid=master_pid, wanted_child=wanted_child, user_memory=user_memory,
                                   shell=shell, join_command=join_command)
    else:
        return run_command_simple(command=command_list, get_output=get_output, shell=shell,
                                  join_command=join_command)


def yield_file(list_of_files):
    # infinite generator
    c = 0
    while True:
        if c == len(list_of_files): c = 0
        yield list_of_files[c]
        c += 1


def compile_cython():
    for f in os.listdir(CYTHON_FOLDER):
        if 'get_non_overlapping_hits.c' in f:
            remove_file(CYTHON_FOLDER + f)
    run_command(f'python {CYTHON_FOLDER}setup_get_non_overlapping_hits.py build_ext --build-lib {CYTHON_FOLDER}')


def cython_compiled():
    if not os.path.exists(CYTHON_FOLDER + 'get_non_overlapping_hits.c'): return False
    return True




def get_combination_ranges(ranges):
    res = set()
    for r in ranges:
        res.update(range(r[0], r[1] + 1))
    return len(res)


def min_max_scale(X, minX, maxX):
    if minX == maxX:
        return 1
    return (X - minX) / (maxX - minX)


def get_hmm_chunk_size(total_profiles, current_chunk_size, max_chunks):
    n_chunks = total_profiles / current_chunk_size
    if n_chunks < max_chunks: return current_chunk_size
    target_chunk_size = round_to_digit(total_profiles / max_chunks)
    return target_chunk_size


def save_metrics(pickle_path, to_pickle):
    with open(pickle_path, 'wb') as handle:
        pickle_dump(to_pickle, handle)


def load_metrics(pickle_path):
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as handle:
            pickled_results = pickle_load(handle)
            return pickled_results


def recalculate_coordinates(env_from, env_to, overlap_value):
    if env_from < env_to + 1:
        start = env_from
        end = env_to + 1
    else:
        start = env_to + 1
        end = env_from

    hit_range = end - start
    hit_overlap = ceil(overlap_value * hit_range / 2)
    start = start + hit_overlap
    end = end - hit_overlap
    return start, end


def create_translation_table(protein_str, base1, base2, base3):
    res = {}
    for i in range(len(protein_str)):
        codon = base1[i] + base2[i] + base3[i]
        res[codon] = protein_str[i]
    return res


def parse_translation_tables(ncbi_tables):
    # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    # ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
    res = {}
    start_recording = False
    base3 = None
    with open(ncbi_tables) as file:
        line = file.readline()
        while line:
            line = line.strip('\n')
            line = line.strip()
            if start_recording:
                if line.startswith('name "') and not re.search('name "SGC', line): table_name = line.split('"')[1]
                if line.startswith('id '):         table_id = line.split()[1].strip(',').strip('"').strip()
                if line.startswith('ncbieaa '):    protein_str = line.split()[1].strip(',').strip('"').strip()
                if line.startswith('-- Base1 '):   base1 = line.split()[2].strip(',').strip('"').strip()
                if line.startswith('-- Base2 '):   base2 = line.split()[2].strip(',').strip('"').strip()
                if line.startswith('-- Base3 '):   base3 = line.split()[2].strip(',').strip('"').strip()
                if base3:
                    translation_table = create_translation_table(protein_str, base1, base2, base3)
                    res[int(table_id)] = {'name': table_name, 'table': translation_table}
                    base3 = None
            if 'Genetic-code-table' in line: start_recording = True
            line = file.readline()
    return res


def translate_protein_seq(sequence, translation_table):
    step = 3
    res = []
    for i in range(0, len(sequence), step):
        codon = sequence[i:i + step]
        try:
            if len(codon) == step:
                aminoacid = translation_table['table'][codon]
                res.append(aminoacid)
        except:
            raise Exception
    return ''.join(res)


def write_translated_fasta(original_fasta_path, translated_fasta_path, translation_table, sample_type):
    with open(translated_fasta_path, 'w+') as file:
        for query, seq in read_protein_fasta_generator(original_fasta_path):
            if sample_type == 'rna':
                seq = seq.replace('U', 'T')
            elif sample_type == 'protein':
                raise Exception
            protein_seq = translate_protein_seq(seq, translation_table)
            chunks = [protein_seq[x:x + 60] for x in range(0, len(protein_seq), 60)]
            chunk_str = '\n'.join(i for i in chunks)
            file.write('>' + query + '\n' + chunk_str + '\n')


def check_sample_type(sample_path):
    res = {}
    c = 0
    with open(sample_path) as file:
        line = file.readline()
        while line and c <= 1000:
            if not line.startswith('>'):
                line = line.strip('\n')
                line = line.strip()
                for residue in line:
                    residue = residue.upper()
                    if residue not in res: res[residue] = 0
                    res[residue] += 1
                c += 1
            line = file.readline()
    sorted_keys = sorted(res.keys())
    if sorted_keys == ['A', 'C', 'G', 'T']: return 'dna'
    if sorted_keys == ['A', 'C', 'G', 'U']:
        return 'rna'
    else:
        return 'protein'


def count_residues(sample_path):
    res = 0
    with open(sample_path) as file:
        line = file.readline()
        while line:
            if not line.startswith('>'):
                line = line.strip('\n')
                line = line.strip()
                res += len(line)
            line = file.readline()
    return res
