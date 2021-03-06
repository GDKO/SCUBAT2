#!/usr/bin/env python

import sys
import re
import argparse
import numpy
import time
import subprocess
from string import maketrans
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio._utils import getattr_str, trim_str
from collections import OrderedDict

#####
# Parameters
#####
DEFAULT_NAME = 'SCUBAT_v2.py'
VERSION = 2.10
DEFAULT_TRANSCRIPT_IDENTITY_CUTOFF = 95
DEFAULT_HSP_IDENTITY_CUTOFF = 95
DEFAULT_EXON_OVERLAP_CUTOFF = 80
DEFAULT_TRANSCRIPT_OVERLAP_CUTOFF = 50
DEFAULT_LIBRARY_INSERT_SIZE = 100
DEFAULT_TRANSCRIPT_COVERAGE_CUTOFF = 90
DEFAULT_TRANSCRIPT_COVERAGE_ONE_CONTIG = 70
DEFAULT_NUMBER_OF_NS = 100

parser = argparse.ArgumentParser(prog=DEFAULT_NAME,description='Scaffold contigs using Transcripts')
parser.add_argument('-b','--blastfile', dest='blastfile', required=True,
                   help='BLAST output in xml format [Required]')
parser.add_argument('-f','--fastafile', dest='fastafile', required=True,
                   help='Contig assembly in fasta format [Required]')
parser.add_argument('-hid','--hsp_identity_cutoff', type=int,
                   help='BLAST HSP identity cutoff [default: %(default)s]',default=DEFAULT_HSP_IDENTITY_CUTOFF)
parser.add_argument('-tid','--transcript_identity_cutoff', type=int,
                   help='BLAST transcript identity cutoff [default: %(default)s]',default=DEFAULT_TRANSCRIPT_IDENTITY_CUTOFF)
parser.add_argument('-tcov','--transcript_coverage_cutoff', type=int,
                   help='Transcript coverage cutoff [default: %(default)s]',default=DEFAULT_TRANSCRIPT_COVERAGE_CUTOFF)
parser.add_argument('-max','--maximum_intron_size', type=int, required=True,
                   help='Maximum intron size [Required]')
parser.add_argument('-eov','--exon_overlap_cutoff', type=int,
                   help='Exon overlap cutoff [default: %(default)s]',default=DEFAULT_EXON_OVERLAP_CUTOFF)                   
parser.add_argument('-tov','--transcript_overlap_cutoff', type=int,
                   help='Transcript overlap cutoff [default: %(default)s]',default=DEFAULT_TRANSCRIPT_OVERLAP_CUTOFF)
parser.add_argument('-lis','--library_insert_size', type=int,
                   help='Library insert size [default: %(default)s]',default=DEFAULT_LIBRARY_INSERT_SIZE)
parser.add_argument('-ns','--number_of_ns', type=int,
                   help='Number of Ns to add when merging [default: %(default)s]',default=DEFAULT_NUMBER_OF_NS)
parser.add_argument('-isr','--intron_size_run', dest='intron_size_run', action='store_true',
                   help='Estimate insert size based on the data')
parser.add_argument('--verbose', dest='verbose', action='store_true',
                   help='Print additional verbose files')
parser.add_argument('--version', action='version', version='%(prog)s v' + str(VERSION))

parser.set_defaults(verbose=False)
parser.set_defaults(intron_size_run=False)

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

#####
# Classes
#####

####
# Class to contain hsps from different hits in one object
####
class Transcript:
    def __init__(self, query_id='<unknown id>',query_length='<unknown length>'):
        self.query_id = query_id
        self.query_length = query_length
        self.hsps = []
        self._items = []
    
    def add(self,hsp):
        self.hsps.append(hsp)
        
    def __iter__(self):
        return iter(self.hsps)
          
    def sort(self, key=None, reverse=False, in_place=True):
        if in_place:
            self.hsps.sort(key=key, reverse=reverse)
        else:
            hsps = self.hsps[:]
            hsps.sort(key=key, reverse=reverse)
            obj = self.__class__(hsps)
            self._transfer_attrs(obj)
            return obj
            
    #print output of the object    
    def __str__(self):
        lines =[]
        qid_line = 'Query: %s' % self.query_id
        qlength_line = 'Length: %i' % self.query_length
        lines.append(qid_line)
        lines.append(qlength_line)
        if not self.hsps:
            lines.append(' HSPs: ?')
        else:
            lines.append(' HSPs: %s  %s  %s  %s  %s  %s  %s  %s  %s  %s' %
                    ('-'*4, '-'*8, '-'*9, '-'*6, '-'*6, '-'*15, '-'*21, '-'*58, '-'*15, '-'*7))
            pattern = '%11s  %8s  %9s  %6s  %6s  %15s  %21s  %58s  %15s  %7s'
            lines.append(pattern % ('#', 'E-value', 'Bit score',  'ID (%)','Span',
                    'Query range', 'Hit range', 'Contig', 'Contig Length', 'Strand'))
            lines.append(pattern % ('-'*4, '-'*8, '-'*9, '-'*6, '-'*6, '-'*15, '-'*21, '-'*58, '-'*15, '-'*7))
            for idx, hsp in enumerate(self.hsps):
                # evalue
                evalue = getattr_str(hsp, 'evalue', fmt='%.2g')
                # bitscore
                bitscore = getattr_str(hsp, 'bitscore', fmt='%.2f')
                # identities
                identity = getattr_str(hsp, 'identity', fmt='%4.2f')
                # alignment length
                aln_span = getattr_str(hsp, 'aln_span')
                # query region
                query_start = getattr_str(hsp, 'query_start')
                query_end = getattr_str(hsp, 'query_end')
                query_range = '[%s:%s]' % (query_start, query_end)
                # max column length is 18
                query_range = trim_str(query_range, 15, '~]')
                # hit region
                hit_start = getattr_str(hsp, 'hit_start')
                hit_end = getattr_str(hsp, 'hit_end')
                hit_range = '[%s:%s]' % (hit_start, hit_end)
                hit_range = trim_str(hit_range, 21, '~]')
                # contig id
                contig_id = getattr_str(hsp, 'hit_id')
                # contig length
                contig_length = getattr_str(hsp, 'hit_length')
                # strand
                strand = getattr_str(hsp, 'strand')
                # append the hsp row
                lines.append(pattern % (str(idx), evalue, bitscore, identity, aln_span,
                        query_range, hit_range, contig_id, contig_length, strand))        
        return '\n'.join(lines)        

####
# Class that contains the information for each hsp
# Named HSPG to avoid conflict with class HSP of biopython module
####
class Hspg:
    def __init__(self,evalue,bitscore,identity,aln_span,query_start,query_end,hit_start,hit_end,hit_id,hit_length, strand):
        self.evalue = evalue
        self.bitscore = bitscore
        self.identity = identity
        self.aln_span = aln_span
        self.query_start = query_start
        self.query_end = query_end
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.hit_id = hit_id
        self.hit_length = hit_length
        self.strand = strand
    def _query_range_get(self):
        return (self.query_start, self.query_end)    
    query_range = property(fget=_query_range_get,
        doc="""Tuple of query start and end coordinates""")

       
#####
# functions
#####

####
# remove hsps if they overlap
# choose the one with the highest bitscore
# (primarly to remove hsps in a pseudogene area in the same contig)
####
def removeOverlappingHsps( hit ):
    sort_key = lambda hsp: hsp.query_range
    hit.sort(key=sort_key, reverse=False, in_place=True)
    k = 0
    while k < len(hit)-1:
        x = range(hit.hsps[k].query_start,hit.hsps[k].query_end)
        y = range(hit.hsps[k+1].query_start,hit.hsps[k+1].query_end)
        d = range(max(x[0], y[0])-1, min(x[-1], y[-1])+1)
        minimum_length = min(len(x),len(y))
        diff = float(len(d)) / (minimum_length+1)    
        if diff > round(float(args.exon_overlap_cutoff)/100,2):
            if hit.hsps[k].bitscore <= hit.hsps[k+1].bitscore:
                hit.hsps.pop(k)
                k -= 1
            else:
                hit.hsps.pop(k+1) 
                k -= 1
        k += 1
	

####
# remove hits if they overlap
# choose the longest one
# (primarly to remove hits in a pseudogene/domain area in other contigs)
####
def removeOverlappingHits (qresult):
    sort_key = lambda hit: len(hit)
    qresult.sort(key=sort_key, reverse=True, in_place=True)
    k = 0
    remove_k = 0
    while k < len(qresult) -1:
        x = range(qresult.hits[k].hsps[0].query_start,qresult.hits[k].hsps[-1].query_end)
        i = k + 1
        done = 0
        while i < len(qresult) and done==0:
            y = range(qresult.hits[i].hsps[0].query_start,qresult.hits[i].hsps[-1].query_end)
            d = range(max(x[0], y[0])-1, min(x[-1], y[-1])+1)
            minimum_length = min(len(x),len(y))
            diff = float(len(d)) / (minimum_length+1)
            if diff > round(float(args.transcript_overlap_cutoff)/100,2):
                if len(x) < len(y):
                    remove_k = 1
                    done = 1
                else:
                    qresult.pop(i)
                    i -= 1
            i += 1
        if remove_k == 1:
            qresult.pop(k)
            k -= 1
            remove_k = 0
        k += 1        

        
        

####
# percentage status calculator
####     
def disp_status(prefix,processed,total):
    if processed and total:
        percent = round(100 * (float(processed)/float(total)),2)
        sys.stdout.write("\r")
        sys.stdout.write(prefix+"%4.2f%% " % percent)
        sys.stdout.flush()


####
# calculate intron size
####
def propable_intron_size(hit_start_first,hit_end_first,orientation_first,contig_length_first,hit_start_second,hit_end_second,orientation_second,contig_length_second):
    if (orientation_first == "+" and orientation_second == "+"):
        p_intron_size = (contig_length_first - hit_end_first) + hit_start_second 
    elif (orientation_first == "+" and orientation_second == "-"):
        p_intron_size = (contig_length_first - hit_end_first) + (contig_length_second - hit_end_second)  
    elif (orientation_first == "-" and orientation_second == "-"):
        p_intron_size = hit_start_first + (contig_length_second - hit_end_second)    
    elif (orientation_first == "-" and orientation_second == "+"):
        p_intron_size = hit_start_first + hit_start_second
    return p_intron_size    

####
# create the path
####
def getPath(nodes, start):
    path = [start]
    while nodes.has_key(path[-1]):
        path.append(nodes[path[-1]])

    return path

####
# return the correct orientation of the contig
####
def getContig(ctgs, node):
    if node[-1] == 'f':
        return ctgs[node[:-2]]
    else:
        return str(Seq(ctgs[node[:-2]]).reverse_complement())

    
#####
# Begin
#####
localtime = time.asctime( time.localtime(time.time()) )
print 'Program started: ' + localtime
start = time.time()
number_of_hits =subprocess.check_output(['grep','-c', '<Iteration>', args.blastfile])


h = open('intron_size.txt', 'w')
print >>h, '#intron_size \t transcript_id'
#debug = open('debug.txt', 'w')

if not args.intron_size_run and args.verbose:
    e = open('filtered_blast_output.txt', 'w')
    v = open('Cyto_table_verbose.sif', 'w')
    p = open('prob_intron.txt', 'w')
    c = open('connection_removed.txt', 'w')
    tr = open('toremove.dict.txt', 'w')
    qcovfile = open('transcript_status.txt', 'w')
    print >>qcovfile, '#transcript_id' + '\t' + 'transcript_coverage' + '\t' + 'identity_percentage' + '\t' + 'status'
    all_connections = open('all_connections.txt', 'w')
    unmapped_transcripts_file = open('unmapped_transcripts.txt', 'w')

if not args.intron_size_run:
    scaffold_file = open('SCUBAT_scaffolds.fasta', 'w')
    confiledict = open('connections.dict.txt', 'w')
    s = open('SCUBAT_stats.txt', 'w')

####
#Process BLAST xml
####
no_hit = []
uninformative = []
informative = []
same_contig = []
n = 0
contig_check_inner = {};
contig_check_outer = {};



for qresult in SearchIO.parse(args.blastfile, 'blast-xml'):
    for hit in qresult:
        removeOverlappingHsps( hit )
    removeOverlappingHits(qresult)
    if len(qresult.hits) == 0:
        no_hit.append(qresult)
        if not args.intron_size_run and args.verbose:
            print >>unmapped_transcripts_file, qresult.id
    elif len(qresult.hits) == 1:
        if len(hit.hsps) == 1:
            uninformative.append(qresult)
        else:    
            same_contig.append(qresult)

        if (hit.hsps[0].hit_start < hit.hsps[-1].hit_end):
            contig_check_start = hit.hsps[0].hit_start
            contig_check_end   = hit.hsps[-1].hit_end
        else:
            contig_check_start = hit.hsps[-1].hit_start
            contig_check_end   = hit.hsps[0].hit_end
        if hit.id in contig_check_inner:
            if contig_check_start < contig_check_inner[hit.id]:
                contig_check_inner[hit.id] = contig_check_end
            elif contig_check_end > contig_check_outer[hit.id]:
                contig_check_outer[hit.id] = contig_check_start
            else:
                pass
        else:               
            contig_check_inner[hit.id] = contig_check_end
            contig_check_outer[hit.id] = contig_check_start
    else:
        transcript = Transcript(qresult.id,qresult.seq_len)
        identity_number = 0
        alignment_hsp = 0
        for hit in qresult:
            for hsp in hit:
                identity_hsp = round(100*(float(hsp.ident_num)/hsp.aln_span),2)
                if identity_hsp >= args.hsp_identity_cutoff:
                    identity_number += hsp.ident_num
                    alignment_hsp += hsp.aln_span
                    hsp_t = Hspg(hsp.evalue,hsp.bitscore,identity_hsp,hsp.aln_span,hsp.query_start,hsp.query_end,hsp.hit_start,hsp.hit_end,hit.id,hit.seq_len,hsp.hit_strand)
                    transcript.add(hsp_t)
        sort_key = lambda hsp: hsp.query_range
        transcript.sort(key=sort_key, reverse=False, in_place=True)
        k = 0
        qcov = []
        while k < len(transcript.hsps):
            qcov.extend(range(transcript.hsps[k].query_start,transcript.hsps[k].query_end+1))
            k += 1
        qcov = list(set(qcov))
        pers = round(100*(float(len(qcov))/transcript.query_length),2)
        # Change to 1 in case of no hsps passing the identity filter
        if alignment_hsp == 0:
            alignment_hsp = 1
        identity = round(100*(float(identity_number)/alignment_hsp),2)
        if (pers >= args.transcript_coverage_cutoff and identity >= args.transcript_identity_cutoff):
            informative.append(transcript)
            if not args.intron_size_run and args.verbose:
                print >>qcovfile, transcript.query_id + '\t' + str(pers) + '\t' + str(identity)
        else:
            if not args.intron_size_run and args.verbose:
                print >>qcovfile, transcript.query_id + '\t' + str(pers) + '\t' + str(identity) + '\t' + 'REMOVED'
    n += 1
    disp_status("Processing BLAST xml...",n,number_of_hits)
print ''

if not args.intron_size_run and args.verbose:
    qcovfile.close()



#####
# Calculate stats for intron sizes and exon overlaps
#####
exon_overlaps = []
intron_sizes = []
n = 0

for qresult in same_contig:
    for hit in qresult:    
        k = 0
        while k < len(hit)-1:
            if (hit.hsps[k+1].query_start-hit.hsps[k].query_end)<1:
                overlap_diff = abs(hit.hsps[k+1].query_start-hit.hsps[k].query_end)+1
                exon_overlaps.append(overlap_diff)
            intron_size = abs(hit.hsps[k+1].hit_start-hit.hsps[k].hit_end)+1
            intron_sizes.append(intron_size)
            h.write(str(intron_size)+ '\t' + hit.query_id +'\n')
            k += 1
    n += 1
    disp_status("Calculating stats...",n,len(same_contig))    
print ''


#h.close()


#####
# Create connections
#####
n = 0
orientation = {}
conn_removed = 0
connections=[]
connections_dict={}
to_remove=[]


for transcript in informative:
    k=0
    if not args.intron_size_run and args.verbose:
        print >>e, transcript
    while k < len(transcript.hsps) -1:
        flag=1
        if (transcript.hsps[k].hit_id == transcript.hsps[k+1].hit_id):
            intron_size = abs(transcript.hsps[k+1].hit_start-transcript.hsps[k].hit_end)+1
            intron_sizes.append(intron_size)
            h.write(str(intron_size)+ '\t' + hit.query_id +'\n')
        else:
            if (transcript.hsps[k].strand == 1 or transcript.hsps[k].strand == '+'):
                transcript.hsps[k].strand='+'
                if transcript.hsps[k].hit_id in contig_check_outer and transcript.hsps[k].hit_end < contig_check_outer[transcript.hsps[k].hit_id]:
                    flag=0
                else:
                    pass
            else:
                transcript.hsps[k].strand='-'
                if transcript.hsps[k].hit_id in contig_check_inner and transcript.hsps[k].hit_start > contig_check_inner[transcript.hsps[k].hit_id]:
                    flag=0
                else:
                    pass                
            if (transcript.hsps[k+1].strand == 1 or transcript.hsps[k+1].strand == '+'):
                transcript.hsps[k+1].strand='+'
                if transcript.hsps[k+1].hit_id in contig_check_inner and transcript.hsps[k+1].hit_start > contig_check_inner[transcript.hsps[k+1].hit_id]:
                    flag=0
                else:
                    pass                
            else:
                transcript.hsps[k+1].strand='-'
                if transcript.hsps[k+1].hit_id in contig_check_outer and transcript.hsps[k+1].hit_end < contig_check_outer[transcript.hsps[k+1].hit_id]:
                    flag=0
                else:
                    pass                
            if not args.intron_size_run:    
                p_intron = propable_intron_size(transcript.hsps[k].hit_start,transcript.hsps[k].hit_end,transcript.hsps[k].strand,transcript.hsps[k].hit_length,transcript.hsps[k+1].hit_start,transcript.hsps[k+1].hit_end,transcript.hsps[k+1].strand,transcript.hsps[k+1].hit_length)
                if args.verbose:
                    p.write(str(p_intron) + '\t' + transcript.query_id + '\n')
                if ((p_intron+args.library_insert_size) < args.maximum_intron_size and flag):
                    if args.verbose:
                        v.write(transcript.hsps[k].hit_id + '\t' + transcript.query_id + '\t' + transcript.hsps[k+1].hit_id + '\t' + transcript.hsps[k].strand + transcript.hsps[k+1].strand+ '\n')
                    if (transcript.hsps[k].strand == "+" and transcript.hsps[k+1].strand == "+"):
                        connections.append(transcript.hsps[k].hit_id + '/f' + '\t' + transcript.hsps[k+1].hit_id + '/f')
                        connections.append(transcript.hsps[k+1].hit_id + '/r' + '\t' + transcript.hsps[k].hit_id + '/r')
                    elif (transcript.hsps[k].strand == "+" and transcript.hsps[k+1].strand == "-"):
                        connections.append(transcript.hsps[k].hit_id + '/f' + '\t' + transcript.hsps[k+1].hit_id + '/r')
                        connections.append(transcript.hsps[k+1].hit_id + '/f' + '\t' + transcript.hsps[k].hit_id + '/r')
                    elif (transcript.hsps[k].strand == "-" and transcript.hsps[k+1].strand == "+"):
                        connections.append(transcript.hsps[k].hit_id + '/r' + '\t' + transcript.hsps[k+1].hit_id + '/f')
                        connections.append(transcript.hsps[k+1].hit_id + '/r' + '\t' + transcript.hsps[k].hit_id + '/f')
                    else:
                        connections.append(transcript.hsps[k].hit_id + '/r' + '\t' + transcript.hsps[k+1].hit_id + '/r')
                        connections.append(transcript.hsps[k+1].hit_id + '/f' + '\t' + transcript.hsps[k].hit_id + '/f')                  
                else:
                    if args.verbose:
                        c.write(transcript.hsps[k].hit_id + '\t' + transcript.query_id + '\t' + transcript.hsps[k+1].hit_id + '\t' + str(p_intron) + '\n')
                    conn_removed += 1
        k += 1
    n += 1
    disp_status("Creating connections...",n,len(informative))

print ''


if args.intron_size_run:
    sys.exit()

if args.verbose:
    e.close()
    p.close()
    v.close()
    c.close()
    unmapped_transcripts_file.close()

transcripts_passed = 0
transcripts_failed = 0

for transcript in informative:
    pass_flag = 0 ;
    transcript_span = {}
    for hsp in transcript:
        #transcript_span[hsp.hit_id] += hsp.aln_span
        transcript_span[hsp.hit_id] = transcript_span.get(hsp.hit_id, 0) + hsp.aln_span 
    for hit_id in transcript_span:
        if ((transcript_span[hit_id]/transcript.query_length)*100) >= DEFAULT_TRANSCRIPT_COVERAGE_ONE_CONTIG:
            pass_flag  = 1 ;
    if pass_flag  == 1:
        transcripts_passed += 1;
    else:
        transcripts_failed += 1;




connections = list(set(connections))

if args.verbose:
    for conn in connections:
        print >>all_connections, conn
    all_connections.close()
    
precs=[]
sucs=[]
for conn in connections:
    pre,suc=conn.split('\t')
    precs.append(pre)
    sucs.append(suc)

precs_duplicates= list (set([x for x in precs if precs.count(x) > 1]))
sucs_duplicates= list (set([x for x in sucs if sucs.count(x) > 1]))

if args.verbose:
    print >>tr, "#precs"
    for trm in precs_duplicates:
        print >>tr, trm
    print >>tr, "#sucs"
    for trm in sucs_duplicates:
        print >>tr, trm
    print >>tr, "#connections removed"

nodes = {}
start_nodes = set()
all_nodes = set()


for conn in connections:
    pre,suc = conn.split('\t')
    if (pre in precs_duplicates or suc in sucs_duplicates):
        if args.verbose:
            print >>tr, conn
    else:
        nodes[pre] = suc
        start_nodes.add(pre)
        all_nodes.add(pre[:-2])
        all_nodes.add(suc[:-2])

if args.verbose:
    tr.close()
start_nodes -= set(nodes.values())

paths = []
contigs = {}

for record in SeqIO.parse(args.fastafile, "fasta") :
    contigs[record.id] = str(record.seq)


string_of_ns = "N" * args.number_of_ns


i = 0;
for start_node in start_nodes:
    path=getPath(nodes, start_node)
    if path[0] > path[-1]:
        i += 1    
        print >>confiledict, path
        scaffold_file.write('>Path_' + str(i) +'_' + str(path) + '\n' + string_of_ns.join([getContig(contigs,node) for node in path]) + '\n')


for contig_id in contigs.keys():
    if contig_id not in all_nodes:
        scaffold_file.write('>' + contig_id + '\n' + contigs[contig_id] + '\n')

confiledict.close()
scaffold_file.close()


# END OF SCAFFOLDING


transcripts_passed_number = len(uninformative) + len(same_contig) + transcripts_passed

# print the stats
s.write('Program was called as:' + '\n')
s.write(DEFAULT_NAME + '\n')
s.write(' -b '    + '\t' + args.blastfile                       + '\n')
s.write(' -f '    + '\t' + args.fastafile                       + '\n')
s.write(' -hid '  + '\t' + str(args.hsp_identity_cutoff)        + '\n')
s.write(' -tid '  + '\t' + str(args.transcript_identity_cutoff) + '\n')
s.write(' -tcov ' + '\t' + str(args.transcript_coverage_cutoff) + '\n')
s.write(' -max '  + '\t' + str(args.maximum_intron_size)        + '\n')
s.write(' -eov '  + '\t' + str(args.exon_overlap_cutoff)        + '\n')
s.write(' -tov '  + '\t' + str(args.transcript_overlap_cutoff)  + '\n')
s.write(' -lis '  + '\t' + str(args.library_insert_size)        + '\n')
s.write(' -ns '   + '\t' + str(args.number_of_ns)               + '\n')
s.write('---------------------------------------------' + '\n')        
s.write('Transcripts with no hits:' + str(len(no_hit))  + '\n')
s.write('Transcripts with 1 hsp:' + str(len(uninformative)) + '\n')
s.write('Transcripts hitting 1 contig with multiple hsps:'  + str(len(same_contig)) + '\n')
s.write('Transcripts hitting multiple contigs:'             + str(len(informative)) + '\n')
s.write('Intron_size Mean:'  + str(numpy.array(intron_sizes).mean(axis=0))   + '\t' + 'Median:' + str(numpy.median(numpy.array(intron_sizes), axis=0))  + '\t' + 'SD:' + str(numpy.array(intron_sizes).std(axis=0))  + '\n')
s.write('Exon_Overlap Mean:' + str(numpy.array(exon_overlaps).mean(axis=0))  + '\t' + 'Median:' + str(numpy.median(numpy.array(exon_overlaps), axis=0)) + '\t' + 'SD:' + str(numpy.array(exon_overlaps).std(axis=0)) + '\n')
s.write('Connections removed based on max intron size:' + str(conn_removed)  + '\n')
s.write('---------------------------------------------' + '\n')
s.write('---------------------------------------------' + '\n')
s.write('Transcripts with no hits:' + str(len(no_hit))  + '\n')
s.write('Transcripts more than ' + str(DEFAULT_TRANSCRIPT_COVERAGE_ONE_CONTIG) + ' in one contig:' + str(transcripts_passed_number) + '\n')
s.write('Transcripts less than ' + str(DEFAULT_TRANSCRIPT_COVERAGE_ONE_CONTIG) + ' in one contig:' + str(transcripts_failed)        + '\n')

s.close()

#####
# Ending prints
#####
localtime = time.asctime( time.localtime(time.time()) )
end = time.time()
print 'Program ended: ' + localtime
seconds = round((end - start),3)
minutes = round((seconds/60),3)
hours = round((minutes/60),3)
print 'SCUBATv2 finished in ' + str(seconds) + ' seconds or ' + str(minutes) + ' minutes or ' + str(hours) + ' hours'
