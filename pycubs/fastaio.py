#!/usr/bin/env python
# coding: utf-8

__author__ = "Author: Guisen Chen; Email: thecgs001@foxmail.com; Date: 2024/05/17"
__all__ = ['FastaIO']
__version__ = "v0.01"

import gzip

def _fastaIO(file):
    """
    Description 
    ----------
    Reads the Fasta format file and returns a generator.
    
    Parameters
    ----------
    file: str
        Read a fasta file path and support compressed files ending in ".gz", or accept a handle of "_io.TextIOWrapper" class.
    """
    
    if isinstance(file, str):
        if file.endswith('.gz'):
            handle = gzip.open(file, 'rt')
        else:
            handle = open(file, 'r')
            
    elif type(file) == '_io.TextIOWrapper':
        handle = file
    else:
        pass
    
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            ID = title.split()[0]
            break
    
    lines = []
    for line in handle:
        if line[0] == ">":
            #yield title, "".join(lines).replace(" ", "").replace("\r", "")
            yield ID, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            ID = title.split()[0]
            continue
        lines.append(line.rstrip())
    #yield title, "".join(lines).replace(" ", "").replace("\r", "")
    yield ID, "".join(lines).replace(" ", "").replace("\r", "")
    handle.close()
    
def fastaIO(file, min_len=0):
    """
    Description 
    ----------
    Reads the Fasta format file and returns a generator.
    
    Parameters
    ----------
    file: str
        Read a fasta file path and support compressed files ending in ".gz", or accept a handle of "_io.TextIOWrapper" class.
        
    min_len: int
        Sequences smaller than this length will be filtered.
    """
    
    if min_len == 0:
        for record in _fastaIO(file):
            yield record
    else:
        for record in _fastaIO(file):
            if len(record[1]) >= min_len: 
                yield record
