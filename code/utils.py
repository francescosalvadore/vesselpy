#!/usr/bin/env python3

import os
import subprocess as sp
import sys
import numpy as np
import json
import math

def verify_openfoam():
    if "FOAM_APPBIN" not in os.environ:
        print('Error! OpenFOAM environment not properly set. Exiting...')
        sys.exit()

def myrun(cmd, shell=True):
    result = sp.run(cmd, shell=shell, capture_output=True, text=True)
    return result.stdout, result.stderr

def print_log(*args, filename="progress.log", end="newline"):
    if end == "newline":
        nl = "\n"
    else:
        nl = ""
    if os.path.exists(filename):
        with open(filename, "a") as f:
            f.write(str(args)+nl)
    else:
        with open(filename, "w") as f:
            f.write(str(args)+nl)

def myround(a, digits=3):
    if isinstance(a, list):
        if len(a) == 0:
            return a
        elif isinstance(a[0], float):
            return [round(x,digits) for x in a]
        elif isinstance(a[0], list):
            a_round = []
            for elem in a:
                elem_rounded =  [round(x,digits) for x in elem]
                a_round.append(elem_rounded)
            return a_round
        else:
            return a
    elif isinstance(a,float):
        return round(a,digits)
    else:
        return a

def decode_range_input(input_string, special="auto"):

    ret = []
    err = ""

    if input_string.strip() == special:
        return special

    if "," in input_string:
        chunks = input_string.split(",")
    else:
        chunks = [input_string]

    for chunk in chunks:
        if ":" in chunk:
            chunk_parts = chunk.split(":")
            if len(chunk_parts) == 3:
                chunk_start = float(chunk_parts[0])
                chunk_end   = float(chunk_parts[1])
                chunk_step  = float(chunk_parts[2])
                max_counter = 1000
                for i in range(max_counter):
                    elem = chunk_start + i*chunk_step
                    if elem > chunk_end:
                        break
                    ret.append(elem)
            else:
                err = f"Chunk {chunk} requires three numbers separated by :"

            pass
        else:
            ret.append(float(chunk))

    return ret #, err

if __name__ == "__main__":
    input_test = "1.,2,5,2.3,4.54"
    input_test = "1.:2.9:0.4"
    input_test = "1.:3:0.4,4.2,7:11:1.2"
    input_test = "5.345"
    input_test = "0:0:10"
    input_test = "-20"
    input_test = "20"
    print(decode_range_input(input_test))
