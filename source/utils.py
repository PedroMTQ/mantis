import requests
from sys import platform,argv
import subprocess
import re
import psutil
import os
import gzip
import sys
from datetime import datetime
from time import sleep,time
from math import ceil,log10
if platform.startswith('win'):    splitter = '\\'
else:                                splitter = '/'
from pathlib import Path
import shutil
import urllib.request as request
from contextlib import closing
from functools import wraps
import csv
import sqlite3
from multiprocessing import Process, current_process,cpu_count, Manager,Pool
from multiprocessing.dummy import Pool as IO_Pool
from zipfile import ZipFile

manager=Manager()
global_queue=manager.list()

def get_slurm_value(wanted_val,regex_pattern):
    res=None
    slurm_job_id = os.environ.get('SLURM_JOBID')
    if slurm_job_id:
        process=subprocess.run('sacct -j '+str(slurm_job_id)+' -o '+wanted_val,shell=True, stdout=subprocess.PIPE)
        wanted=re.search(regex_pattern,str(process.stdout))
        if wanted: res=wanted.group()
    return res


def check_environment_cores():
    environment_cores=get_slurm_value('AllocCPUS',re.compile('\d+'))
    if environment_cores:
        if int(environment_cores):
            #print('Cores allocated by slurm:', environment_cores)
            return int(environment_cores)
        else:
            environment_cores=cpu_count()
            #print('Cores allocated:', environment_cores)
            return int(environment_cores)
    else:
        environment_cores = cpu_count()
        #print('Cores allocated:', environment_cores)
        return int(environment_cores)

def check_process_overhead():
    process = psutil.Process()
    ram_consumption_gb = process.memory_info().rss / (1024 * 1024 * 1024)
    del process
    return ram_consumption_gb

def check_available_ram():
    required_ram=get_slurm_value('ReqMem',re.compile('\d+G[nc]'))
    if required_ram:
        if 'Gc' in  required_ram:
            required_ram = required_ram.replace('Gc','')
            ram_per_core=int(required_ram)
            available_ram = environment_cores * ram_per_core
        if 'Gn' in required_ram:
            required_ram = required_ram.replace('Gn','')
            available_ram = int(required_ram)
        #print('RAM allocated by slurm:',available_ram,'GB')
    else:
        #im not sure if this is the best implementation
        available_ram= psutil.virtual_memory().available/(1024*1024*1024)
        #print('RAM allocated:',available_ram,'GB')
    return available_ram

#this is the ovehead generated per each python process. Multiprocessing generates new python processes per Process
process_overhead = check_process_overhead()
environment_cores = check_environment_cores()
available_ram = check_available_ram()

worker_per_core = 1

mantis_folder = os.path.abspath(os.path.dirname(__file__)).split(splitter)[0:-1]
mantis_folder = splitter.join(mantis_folder)+splitter


def estimate_number_workers_annotation(n_chunks=0,
                                       n_hmms=0,
                                       default_workers=None,
                                       user_cores=None,
                                       minimum_jobs_per_worker=1,
                                       percentage_allowed_overhead=0.9,
                                       split_sample=False
                                       ):
    if default_workers: return default_workers
    if user_cores:
        environment_workers=user_cores*worker_per_core
    else:
        environment_workers=environment_cores*worker_per_core
    #allowed number of workers in order to preserve the maximum allowed overhead reserved for process spawning
    maximum_overhead_workers=ceil((available_ram*percentage_allowed_overhead)/process_overhead)
    #maximum amount of workers that could be spawning considering the ram required for the annotation and the overhead associated with spawning a process
    #total times we need to run hmmer
    total_hmmer_runs=n_chunks*n_hmms
    if split_sample: total_hmmer_runs=environment_workers
    return ceil(min([total_hmmer_runs,maximum_overhead_workers,environment_workers])/minimum_jobs_per_worker)





def estimate_chunk_size(total_n_seqs,
                        annotation_workers,
                        chunk_size=None,
                        minimum_chunk_size=200,
                        ):
    '''
    this is a double edged sword, splitting into too many small chunks generates too much overhead
    not splitting into enough chunks leads to idle processes
    '''
    if chunk_size: return chunk_size
    #maximize chunk size taken into account the amount of workers and number of seqs
    potential_chunk_size= round_to_digit(total_n_seqs/annotation_workers)
    if available_ram<=4: maximum_chunk_size=500
    elif available_ram>4 and available_ram<=10: maximum_chunk_size=1000
    elif available_ram>10 and available_ram<=30: maximum_chunk_size=2000
    elif available_ram>30 and available_ram<=50: maximum_chunk_size=5000
    elif available_ram>50 and available_ram<=70: maximum_chunk_size=10000
    else: maximum_chunk_size=20000
    if potential_chunk_size<minimum_chunk_size:     return minimum_chunk_size
    elif potential_chunk_size>maximum_chunk_size:   return maximum_chunk_size
    else:                                           return potential_chunk_size

def estimate_number_workers_split_sample(minimum_worker_load,
                                         len_protein_seqs,
                                         percentage_allowed_overhead=0.8):
    estimate_chunk_workers = ceil(len_protein_seqs / minimum_worker_load)
    #we only reserve a certain percentage for the overhead associated with spawning processes
    maximum_overhead_workers=int((available_ram*percentage_allowed_overhead)/process_overhead)
    return min([maximum_overhead_workers,estimate_chunk_workers])

def estimate_number_workers_setup_database(n_hmms_to_setup,
                                           percentage_allowed_overhead=0.8,
                                           minimum_jobs_per_worker=3):
    #0.1GB RAM per hmm setup
    ram_per_setup=0.1
    maximum_overhead_workers=ceil((available_ram*percentage_allowed_overhead)/process_overhead)
    maximum_workers= ceil(available_ram/(ram_per_setup+process_overhead))
    return ceil(min([n_hmms_to_setup,maximum_workers,maximum_overhead_workers])/minimum_jobs_per_worker)

def estimate_number_workers_process_output(n_chunks,domtblout_per_chunks=1):
    return min([environment_cores*worker_per_core,n_chunks*domtblout_per_chunks])




def is_ec(enz_id,required_level=3):
    if enz_id:
        ec_pattern = re.compile('^\d+\.\d+\.\d+(\.(-|\d+|([a-zA-Z]\d+)))?')
        if re.search(ec_pattern, enz_id):
            enz_id_copy=str(enz_id).replace('.-','')
            if len(enz_id_copy.split('.'))>=required_level:
                return True
    return False


def find_ecs_legacy(string_to_search,required_level=3):
    res=set()
    ec_pattern = re.compile('\d+\.\d+\.\d+(\.(-|\d+|([a-zA-Z]\d+)))?')
    search= re.finditer(ec_pattern,string_to_search)
    for i in search:
        ec=i.group()
        if string_to_search[i.span()[0]-1]!='.':
            if is_ec(ec,required_level):
                if ec[-2:]=='.-': ec=ec[0:-2]
                res.add(ec)
    return res

def find_ecs(string_to_search, required_level=3):
    res = set()
    # greedy match of confounders
    ec_pattern = re.compile('\d(\.(-|\d{1,3}|([a-zA-Z]\d{1,3}))){2,3}')
    search = re.finditer(ec_pattern, string_to_search)
    for i in search:
        ec = i.group()
        passed = False
        start = i.span()[0]
        end = i.span()[1]
        if len(string_to_search) > end:
            if string_to_search[start - 1] != '.' and\
                    string_to_search[end] != '.' \
                    and not re.match('\.|[a-zA-Z]|\d{1,3}', string_to_search[end]) and not re.match('-', string_to_search[end]):
                passed = True
        else:
            if string_to_search[start - 1] != '.':
                passed = True
        if passed:
            if ec.count('.') >= required_level - 1:
                if ec.count('.') + 1 - ec.count('-') >= required_level:
                    res.add(ec)
    return res

def find_tcdb(string_to_search):
    res=set()
    tc_pattern = re.compile('\(TC\s\d\.[A-Z\-](\.(\d+|\-)){1,2}\)')
    search= re.finditer(tc_pattern,string_to_search)
    for i in search:
        res.add(i.group())
    return res

def find_ko(string_to_search):
    res=set()
    #I could do upper and lower case but since it's only one letter, it's not very safe...
    pattern = re.compile('K\d{4,}')
    search= re.finditer(pattern,string_to_search)
    for i in search:
        res.add(i.group())
    return res

def find_pfam(string_to_search):
    res=set()
    duf = re.compile('(DUF|duf)\d+')
    pfam = re.compile('((U|u)?PF|pf)\d+')
    for pattern in [duf,pfam]:
        search= re.finditer(pattern,string_to_search)
        for i in search:
            res.add(i.group())
    return res

def find_cog(string_to_search):
    res=set()
    pattern = re.compile('(COG|cog)\d+')
    search= re.finditer(pattern,string_to_search)
    for i in search:
        res.add(i.group())
    return res

def find_go(string_to_search):
    res=set()
    pattern = re.compile('GO\d+')
    search= re.finditer(pattern,string_to_search)
    for i in search:
        res.add(i.group())
    return res


def convert_time(user_hours):
    if not user_hours:
        return None
    elif user_hours=='None':
        return None
    elif isinstance(user_hours,int) or isinstance(user_hours,float):
        return ceil(float(user_hours))
    else:
        return user_hours

def get_sample_metrics(target_sample):
    total_sample_length=0
    seqs_length=[]
    with open(target_sample) as file:
        line=file.readline()
        #skip first line with >
        line=file.readline()
        seq_c =0
        while line:
            if line[0]!='>':
                line=line.strip('\n')
                total_sample_length+=len(line)
                seq_c+=len(line)
            else:
                seqs_length.append(seq_c)
                seq_c = 0
            line=file.readline()
        seqs_length.append(seq_c)
    return total_sample_length,total_sample_length/len(seqs_length),len(seqs_length)

def get_seqs_count(target_sample):
    total_seqs=0
    with open(target_sample) as file:
        line=file.readline()
        while line:
            if line[0]=='>': total_seqs+=1
            line=file.readline()
    return total_seqs

def get_folder(file_path):
    return splitter.join(file_path.split(splitter)[0:-1]) + splitter

def read_profile(hmm_file):
    res=[]
    line=hmm_file.readline()
    while line:
        #print(repr(line),line[0:2])
        res.append(line)
        if '//' in line [0:2]:
            return res
        line=hmm_file.readline()

def get_hmm_in_folder(folder_path):
    if os.path.exists(folder_path):
        if folder_path[0:2] != 'NA':
            files =os.listdir(folder_path)
            for f in files:
                if '.hmm' in f[-4:]:
                    return folder_path+f
    return None



def get_hmm_profile_count(hmm_path):
    res=0
    with open(hmm_path,'r') as hmm_file:
        line=hmm_file.readline()
        while line:
            if 'NAME' in line: res+=1
            line=hmm_file.readline()
    return res

def compile_hmm_chunks_path(hmm_path):
    #this will check if the hmm have been split into chunks, if so it will retrieve the chunks instead
    res=[]
    hmm_folder=get_folder(hmm_path)
    if 'chunks' in os.listdir(hmm_folder):
        for chunk_file in os.listdir(add_slash(hmm_folder)+'chunks/'):
            if chunk_file[-4:]=='.hmm': res.append(add_slash(hmm_folder)+'chunks/'+chunk_file)
        return res
    else:
        return [hmm_path]


def round_to_digit(number):
    number_str=str(int(number))
    res=number_str[0:1]
    to_add=1
    if len(number_str)>1:
        if int(number_str[1])>=5:
            to_add+=1
            res+='5'
    for i in range(len(number_str)-to_add):
        res+='0'
    return int(res)


def chunk_generator(to_chunk, chunk_size):
    for i in range(0, len(to_chunk), chunk_size):
        yield to_chunk[i:i + chunk_size]

def chunk_generator_load_balanced(list_ordered,chunk_size,time_limit=60):
    '''
    :param list_chunks: list of keys which will correpond to sequences IDs . list should have been ordered by sequence length
    :param time_limit: when dealing with a lot of sequences, this method may take a while, if it takes too long we raise TimeoutError
    #1M sequences take up about 80-90 seconds to create the chunks, which is still kinda okay. This should be below most metagenomic samples
    '''
    n_chunks=ceil(len(list_ordered)/chunk_size)
    res=[]
    direction_chunks={}
    for i in range(n_chunks):
        res.append([])
        direction_chunks[i]=True
    chunk_index=0
    start=time()
    while list_ordered:
        if time()-start>time_limit: raise TimeoutError
        if direction_chunks[chunk_index]:
            chunk_val=list_ordered.pop(0)
            direction_chunks[chunk_index]=False
        else:
            chunk_val=list_ordered.pop(-1)
            direction_chunks[chunk_index]=True
        res[chunk_index].append(chunk_val)
        if chunk_index==n_chunks-1: chunk_index=0
        else: chunk_index+=1
    return res


def is_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        line=f.readline()
        if line[0]=='>':
            return True
    return False


def add_slash(path_folder):
    if not path_folder: return path_folder
    if path_folder[-1]!=splitter: return path_folder+splitter
    return path_folder

def get_path_level(path,level=1,remove_extension=False):
    temp=path.strip(splitter)
    res=temp.split(splitter)
    res= res[-level]
    if level==1 and remove_extension: res=res.split('.')[0]
    return res

def print_cyan(*args,flush=False,file=None):
    for i in args:  print('# ', i,flush=flush,file=file)
def red(text,flush=False,file=None):
    print('\033[31m', text, '\033[0m',flush=flush, sep='',file=file)
def green(text,flush=False,file=None):
    print('\033[32m', text, '\033[0m',flush=flush, sep='',file=file)
def yellow(text,flush=False,file=None):
    print('\033[33m', text, '\033[0m',flush=flush, sep='',file=file)

def timeit_function(f):
    @wraps(f)
    def wrapper(self,*args, **kwargs):
        start_time=time()
        res = f(self,*args, **kwargs)
        print('This function',f.__name__,'took',time()-start_time,'seconds to run')
        return res
    return wrapper

def timeit_class(f):
    @wraps(f)
    def wrapper(self,*args, **kwargs):
        res = f(self,*args, **kwargs)
        datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print_cyan('This '+type(self).__name__+' process finished running at ' + datetime_str + ' and took ' + str(int(time() - self.start_time)) + ' seconds to complete.', flush=True, file=self.redirect_verbose)
        return res
    return wrapper

def uncompress_archive(source_filepath, extract_path=None,block_size=65536, remove_source=False, stdout_file=None):
    file_name=source_filepath.split(splitter)[-1]
    dir_path=splitter.join(source_filepath.split(splitter)[0:-1])
    if not extract_path: extract_path=dir_path
    if '.tar' in file_name:
        unpack_archive(source_file=source_filepath, extract_dir=extract_path, remove_source=remove_source, stdout_file=None)
    #only for files
    elif '.gz' in file_name:
        gunzip(source_filepath=source_filepath,dest_filepath=extract_path,block_size=block_size,remove_source=remove_source,stdout_file=stdout_file)
    elif '.zip' in file_name:
        unzip_archive(source_file=source_filepath, extract_dir=extract_path, remove_source=remove_source, stdout_file=None)
    else:
        print('Incorrect format! ',source_filepath,flush=True,file=stdout_file)

#this unzips to the same directory!
def gunzip(source_filepath, dest_filepath=None, block_size=65536,remove_source=False,stdout_file=None):
    print('Gunzipping ',source_filepath,flush=True,file=stdout_file)
    if not dest_filepath:
        dest_filepath=source_filepath.strip('.gz')
    if os.path.isdir(dest_filepath):
        file_name=source_filepath.split(splitter)[-1].replace('.gz','')
        dest_filepath= add_slash(dest_filepath)+file_name
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

def unpack_archive(source_file,extract_dir,remove_source=False,stdout_file=None):
    print('Unpacking ',source_file,flush=True,file=stdout_file)
    shutil.unpack_archive(source_file,extract_dir=extract_dir)
    if remove_source: os.remove(source_file)

def unzip_archive(source_file,extract_dir,remove_source=False,stdout_file=None):
    print('Unzipping ',source_file,flush=True,file=stdout_file)
    with ZipFile(source_file, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)
    if remove_source: os.remove(source_file)


def download_file_http(url,file_path,c):
    if c>5:
        download_file_http_failsafe(url,file_path)
    else:
        with requests.get(url, stream=True) as r:
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)

#slower but safer
def download_file_http_failsafe(url,file_path):
    with requests.Session() as session:
        get = session.get(url, stream=True)
        if get.status_code == 200:
            with open(file_path, 'wb') as f:
                for chunk in get.iter_content(chunk_size=1024):
                    f.write(chunk)


def download_file_ftp(url,file_path):
    with closing(request.urlopen(url)) as r:
        with open(file_path, 'wb') as f:
            shutil.copyfileobj(r, f)

def download_file(url,output_folder='',stdout_file=None,retry_limit=10):
    file_path = output_folder + url.split('/')[-1]
    try:
        target_file = request.urlopen(url)
    except:
        print('Cannot download target url',url)
        return
    target_size= target_file.info()['Content-Length']
    transfer_encoding = target_file.info()['Transfer-Encoding']
    if target_size: target_size=int(target_size)
    if os.path.exists(file_path):
        if transfer_encoding=='chunked':
            return
        elif os.stat(file_path).st_size == target_size:
            print('Not downloading from '+url+' since file was already found!',flush=True,file=stdout_file)
            return
        else: os.remove(file_path)
    print('Downloading from ' + url + '. The file will be kept in '+output_folder,flush=True,file=stdout_file)
    c=0
    while c<=retry_limit:
        if 'ftp' in url:    download_file_ftp(url,file_path)
        else:               download_file_http(url,file_path,c)
        if transfer_encoding=='chunked': return
        if os.stat(file_path).st_size == target_size: return
        c+=1
    print('Did not manage to download the following url correctly:\n'+url)
    raise Exception

def concat_files(output_file,list_file_paths,stdout_file=None):
    print('Concatenating files into ',output_file,flush=True,file=stdout_file)
    with open(output_file,'wb') as wfd:
        for f in list_file_paths:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

def merge_profiles(folder_path,output_file,stdout_file=None):
    print('Merging profiles in ',folder_path,flush=True,file=stdout_file)
    main_folder=get_folder(output_file)
    old_files=os.listdir(main_folder)
    for file in old_files:
        if '.hmm' in file:
            os.remove(main_folder+file)
    list_dir= os.listdir(folder_path)
    profiles = [folder_path+splitter+i for i in list_dir if '.hmm' in i.lower()]
    concat_files(output_file,profiles,stdout_file=stdout_file)
    if os.path.exists(folder_path):      shutil.rmtree(folder_path)

def move_file(source_file,dest_file):
    if not os.path.isdir(dest_file):
        if os.path.exists(dest_file): os.remove(dest_file)
    try:
        os.rename(source_file,dest_file)
    except:
        from shutil import move
        move(source_file,dest_file)

def remove_file(source_file):
    if os.path.exists(source_file): os.remove(source_file)


def run_command_simple(command, get_output=False,stdout_file=None):
    if get_output:
        process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif stdout_file:
        process = subprocess.run(command, stdout=stdout_file, stderr=stdout_file)
    else:
        process = subprocess.run(command)
    return process

def get_available_ram_percentage(worker_status,user_memory=None):
    workers_ram = 0
    worker_ids = list(worker_status.keys())
    for worker_id in worker_ids:
        try:
            workers_ram += psutil.Process(worker_id).memory_info().rss
        except psutil.NoSuchProcess:
            worker_status.pop(worker_id)

    workers_ram/=1024**3
    if user_memory:
        workers_ram_percent=100*workers_ram/user_memory
    else:
        workers_ram_percent=100*workers_ram/available_ram
    return 100-workers_ram_percent

def get_child_workers(master_pid,wanted_name=None):
    '''
    master_pid is the main mantis process (lvl0) which spawns all other multiprocessing processes (lvl1)
    these processes (lvl1) will then run a X command which in turn are themselves processes (lvl2)
    so what we want are not the lvl1 processes but the lvl2 processes. We can know these by their process_name / wanted_name
    '''
    res=set()
    #if not psutil.pid_exists(master_pid): return res
    try:
        master_process=psutil.Process(pid=master_pid)
        children = master_process.children(recursive=True)
        for process in children:
            try:
                if wanted_name:
                    if wanted_name==process.name():
                        res.add(process.pid)
                else:
                    res.add(process.pid)
            except (psutil.NoSuchProcess,AttributeError,FileNotFoundError):
                pass
    except (psutil.NoSuchProcess):
        return res
    return sorted(res)

def count_running_workers(worker_status):
    res=0
    for w in worker_status:
        if worker_status[w]: res+=1
    return res

def get_workers_status(workers_ids):
    res={}
    for worker_pid in workers_ids:
        try:
            process=psutil.Process(pid=worker_pid)
            if process.status()=='running':
                res[worker_pid]=True
            else:
                res[worker_pid]=False
        except psutil.NoSuchProcess: pass
    return res

def kill_workers(master_pid,worker,worker_pid):
    sleep(5)
    try:
        if psutil.pid_exists(master_pid):
            psProcess = psutil.Process(pid=master_pid)
            print('###### Ran out of memory, either increase available memory or try to manually decrease chunk size. Quitting now! ######',flush=True)
            psProcess.kill()
            os.waitpid(master_pid, os.WEXITED)
    except psutil.NoSuchProcess: pass
    try:
        worker.kill()
        os.waitpid(worker_pid, os.WEXITED)
    except psutil.NoSuchProcess: pass

def resume_workers(child_status,n_workers):
    ordered_ids=sorted(child_status.keys())
    if count_running_workers(child_status)==n_workers: return
    for worker_pid in ordered_ids:
        try:
            process = psutil.Process(pid=worker_pid)
            if n_workers>0:
                process.resume()
                child_status[worker_pid]=True
                n_workers-=1
            else:
                process.suspend()
                child_status[worker_pid]=False
        except psutil.NoSuchProcess:
            child_status.pop(worker_pid)


def suspend_workers(child_status,n_workers):
    ordered_ids=sorted(child_status.keys(),reverse=True)
    if not n_workers: return
    if len(child_status)-count_running_workers(child_status)==n_workers: return
    for worker_pid in ordered_ids:
        try:
            process = psutil.Process(pid=worker_pid)
            if n_workers>0:
                process.suspend()
                child_status[worker_pid]=False
                n_workers-=1
            else:
                process.resume()
                child_status[worker_pid]=True
        except psutil.NoSuchProcess:
            child_status.pop(worker_pid)

#in order to not run out of memory when running HMMER we control how the processes are running, suspending them if needed
def run_command_managed(command,master_pid, get_output=False,stdout_file=None,wanted_child=None,user_memory=None):
    if get_output:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif stdout_file:
        process = subprocess.Popen(command, stdout=stdout_file, stderr=stdout_file)
    else:
        process = subprocess.Popen(command)
    process_pid=process.pid
    psProcess = psutil.Process(pid=process_pid)
    while psutil.pid_exists(process_pid):
        if not psutil.pid_exists(master_pid):
            kill_workers(master_pid,process,process_pid)
        if psProcess.status()=='zombie':
            process.communicate()
            return process
        child_workers = get_child_workers(master_pid=master_pid, wanted_name=wanted_child)
        child_status = get_workers_status(child_workers)
        available_ram_percentage=get_available_ram_percentage(child_status,user_memory)
        #if we have more than 20% ram available we let it run
        if available_ram_percentage>=20:
            resume_workers(child_status,n_workers=len(child_status))
        #if we have only 10-20% ram available, we suspend half the processes
        elif available_ram_percentage>=10 and available_ram_percentage<20:
            workers_to_resume=ceil(len(child_status)/2)
            if len(child_status)%2:workers_to_suspend=workers_to_resume-1
            else:workers_to_suspend=workers_to_resume
            resume_workers(child_status,n_workers=workers_to_resume)
            if count_running_workers(child_status)>1:
                suspend_workers(child_status,n_workers=workers_to_suspend)
        #if there's not enough memory to run even a single process, we crash execution
        elif available_ram_percentage < 5:
            kill_workers(master_pid,process,process_pid)
        #if we have less than 10% we suspend all but 1
        else:
            #if we have only 1 running worker we let it run
            if count_running_workers(child_status)<1:
                resume_workers(child_status, n_workers=1)
            #if we have more than one we suspend all but one
            elif count_running_workers(child_status) >1:
                workers_to_suspend = len(child_status)-1
                suspend_workers(child_status,workers_to_suspend)
                resume_workers(child_status, n_workers=1)
        sleep(1)
    return process



def run_command( command, get_output=False,stdout_file=None,master_pid=None,wanted_child=None,user_memory=None):
    command_list=command.split()
    if master_pid:
        return run_command_managed(command=command_list,get_output=get_output,stdout_file=stdout_file,master_pid=master_pid,wanted_child=wanted_child,user_memory=user_memory)
    else:
        return run_command_simple(command=command_list,get_output=get_output,stdout_file=stdout_file)

def yield_file(list_of_files):
    #infine generator
    c=0
    while True:
        if c==len(list_of_files): c=0
        yield list_of_files[c]
        c+=1

def compile_cython():
    cython_folder=mantis_folder+'source'+splitter+'cython_src'+splitter
    for f in os.listdir(cython_folder):
        if 'get_non_overlapping_hits.c' in f:
            remove_file(cython_folder+f)
    run_command('python '+cython_folder+ 'setup_get_non_overlapping_hits.py build_ext '+
                '--build-lib '+cython_folder)

def cython_compiled():
    cython_folder=mantis_folder+'source'+splitter+'cython_src'+splitter
    if not os.path.exists(cython_folder+'get_non_overlapping_hits.c'): return False
    return True

def is_picklable(obj):
    import pickle
    try:
        pickle.dumps(obj)
    except pickle.PicklingError:
        return False
    return True

if __name__ == '__main__':
    compile_cython()