import os
from multiprocessing import cpu_count
from pathlib import Path
import psutil

SERVICE_NAME = 'mantis_pfa'
DEBUG = int(os.environ.get('DEBUG', '0'))
PSUTIL_EXCEPTIONS = (psutil.NoSuchProcess, AttributeError, FileNotFoundError)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STATIC_DATA = os.path.join(ROOT, 'src', 'data')
DATA = os.getenv(os.path.join(ROOT,  'data'), 'MANTIS_DATA')
Path(DATA).mkdir(parents=True, exist_ok=True)

def check_environment_cores():
    res = cpu_count()
    return int(res)


def check_process_overhead():
    process = psutil.Process()
    ram_consumption_gb = process.memory_info().rss / (1024 * 1024 * 1024)
    del process
    return ram_consumption_gb


def check_available_ram():
    res = psutil.virtual_memory().available / (1024 * 1024 * 1024)
    return res



# this is the ovehead generated per each python process. Multiprocessing generates new python processes per Process
PROCESS_OVERHEAD = check_process_overhead()
ENVIRONMENT_CORES = check_environment_cores()
AVAILABLE_RAM = check_available_ram()

# typically should be one worker per physical core, however increasing this further is possible with context switch. Be careful doing so though, as context switching also impacts execution performance.
# having too many workers per core may also demand too much read/writing which can tax the system (which can impact other users of the same system)
WORKER_PER_CORE = int(os.environ.get('MANTIS__WORKER_PER_CORE', '1'))
PERCENTAGE_ALLOWED_OVERHEARD = float(os.environ.get('MANTIS__PERCENTAGE_ALLOWED_OVERHEARD', '0.9'))
MINIMUM_JOBS_PER_WORKER = float(os.environ.get('MANTIS__MINIMUM_JOBS_PER_WORKER', '1'))


DEFAULT_CONFIG = os.path.join(STATIC_DATA, 'config', 'MANTIS.cfg')
CYTHON_FOLDER = os.path.join(ROOT, 'src', 'cython_src')
INSERTION_BATCH_SIZE = 5000
