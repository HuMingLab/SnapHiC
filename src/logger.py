
import sys
import time
import os

class Logger(object):
    def __init__(self, logfile, rank = 0, verbose_threshold = 0):
        self.terminal = sys.stdout
        self.log = open(logfile, "a")
        self.verbose_threshold = verbose_threshold
        self.rank = rank

    def write(self, message, append_time = True, verbose_level = 0, allow_all_ranks = False):
        if append_time:
            message = f'{time.ctime()}: {message}'
        if message[-1] != '\n':
            message += '\n'
        if verbose_level <= self.verbose_threshold:
            if self.rank == 0 or allow_all_ranks:
                self.terminal.write(message)
                self.log.write(message)

    def flush(self):
        self.log.flush()
        self.terminal.flush()

    def set_rank(self, rank):
        self.rank = rank

    def dump_args(self, args):
        self.write("Running snapHiC with following arguments")
        for k, v in vars(args).items():
            self.write(f'\t{k}={v}', append_time = False)
