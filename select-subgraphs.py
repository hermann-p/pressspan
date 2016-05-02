#!/usr/bin/env python3

import re, sys, os, shutil
from subprocess import call

def isValidRegion(region):
    return len(filter(lambda x: len(x) > 0, re.split('[:-]', region))) == 3    

def read_frag(frag):
    chrom, start, end = filter(lambda x: len(x) > 0, re.split('[:-]', frag))
    strandiness = '-' if frag[0] == '-' else '+'
    return {'dir': strandiness,
            'chr': chrom,
            'start': start,
            'end': end}

def is_in_region(region, frag):
    d = region['dir']
    start = region['start']
    end = region['end']
    chrom = region['chr']
    if frag['chr'] == chrom and frag['dir'] == d:
        if frag['start'] >= start and frag['end'] <= end:
            return True
        elif frag['start'] < start and frag['end'] > end:
            return True

def is_in_any_region(regions, el):
    def touples(regs, fs):
        for r in regs:
            for f in fs:
                yield (r,f)
    plotfile, frags, chroms = el.split('\t')
    frags = map(read_frag, frags[1:-1].split(','))
    return reduce(lambda x,y: x or y, [is_in_region(r,f) for r,f in touples(regions, frags)])

def matching_files(filename, regions):
    region_maps = map(read_frag, regions)
    with open(filename, 'r') as inFile:
        matches = filter(lambda x: is_in_any_region(region_maps, x), inFile)
    return map(lambda x: x.split('\t')[0], matches)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: %s region [region region ...] [--multistrand|--circular,--copy]')
        print('\n--multitrand or --circular will only look for respective events')
        print('--copy                       copy findings to subfolder "myfindings/"')
        print('--viz                        create plots with GraphViz')
        exit(1)
    
    regions = filter(lambda x: x[:2] != '--', sys.argv[1:])
    flags = filter(lambda x: x[:2] == '--', sys.argv[1:])

    if '--circular' in flags and '--multistrand' in flags:
        flags = filter(lambda x: x != '--circular' and x != '--multistrand', flags)
    
    if not reduce(lambda x,y: x and y, map(isValidRegion, regions)):
        print('Error: Invalid region definitions.')
        print('Required: [+/-][chromosome-name]:[start-position]-[end-position]')
        print('          -X:100-21313')
        exit(-1)

    results = matching_files('pressspan.log', regions)

    if '--circular' in flags:
        results = filter(lambda x: x.startswith('circ'), results)
    elif '--multistrand' in flags:
        results = filter(lambda x: x.startswith('mult'), results)

    print("Found %d results matching your criteria" % (len(results)))

    if '--copy' in flags and len(results) > 0:
        result_dir = 'myfindings'
        print('Creating folder %s and copying files...' % (result_dir))
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        for rfile in results:
            rfile += '.dot'
            prefix = 'circulars/' if rfile.startswith('circ') else 'multis/'
            try:
                shutil.copyfile(prefix + rfile, result_dir + '/' + rfile)
            except IOError:
                print('ERROR: Could not copy file %s%s' %(prefix, rfile))

    if '--viz' in flags:
        result_dir = 'mygraphs'
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        for rfile in results:
            prefix = 'circulars/' if rfile.startswith('circ') else 'multis/'
            try:
                call(['dot', '-Teps',
                      prefix + rfile + '.dot',
                      '-o', result_dir + '/' + rfile + '.eps'])
            except:
                print('ERROR: Could not complete call to dot. Is GraphViz installed?')
