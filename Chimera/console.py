import os
from datetime import datetime
import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def nominal(message, verbose):
    if verbose is True:
        time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        status = "\r" + bcolors.OKGREEN + ("[~] (" + str(time) + ")  " + str(message)) + bcolors.ENDC
        sys.stdout.write(status)
        sys.stdout.flush()

def error(message, verbose):
    if verbose is True:
        time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        status = bcolors.FAIL + ("\n[X] (" + str(time) + ")  " + str(message)) + bcolors.ENDC
        print(status)

def event(message, verbose):
    if verbose is True:
        time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        status = bcolors.OKBLUE + ("\n[!] (" + str(time) + ")  " + str(message)) + bcolors.ENDC
        print(status)