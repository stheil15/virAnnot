#!/usr/bin/env python

__version_major__ = 1
__version_minor__ = 0
__revision__ = 1
__build__ = ""

import optparse
import sys, os

p = optparse.OptionParser(description = """demultiplex-BLAST-results: split a
XML-formatted BLAST results file into multiple XML-formatted files based on
redirection rules""")

g = optparse.OptionGroup(p, "input and output")

g.add_option("-f", "--fof", dest = "fof", metavar = "FILENAME",
    help = "File with query ID list.")
g.add_option("-i", "--input", dest = "input_fn", metavar = "FILENAME",
    help = "File with query ID list.")

p.add_option_group(g)

p.add_option("-v", "--version", dest = "display_version", action = "store_true", default = False,
    help = "Display the program version and exit")

(p, a) = p.parse_args()

def error (msg):
    print >>sys.stderr, "ERROR: %s" % msg
    sys.exit(1)

if (p.display_version):
    print "%s.%srev%s" % (__version_major__, __version_minor__, __revision__)
    sys.exit(0)

try:
    import amara
    from amara import bindery
except:
    error("The Amara library is required; see http://wiki.xml3k.org/Amara for information.")

if p.fof is None:
    error("At least one list file.")

XPATH = re.compile("(?P<xpath>(?:query|hit|HSP)(?:\.[A-Z][A-Za-z_\-]+)+)")
wrap_xpath = lambda x: XPATH.sub("cast(\g<xpath>)", x)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
print('reading xml...')
print(p.input_fn)
try:
    # note that we do not validate the XML file against
    # http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd
    i = bindery.parse(p.input_fn, standalone = True)

except amara.ReaderError as msg:
    error(str(msg))
print('done.')
n_queries = len(i.BlastOutput.BlastOutput_iterations.Iteration)
print(str(n_queries) + ' queries read.')
# create an envelop from the input XML BLAST
# results containing only the header
envelop = bindery.parse(p.input_fn, standalone = True)
envelop.BlastOutput.BlastOutput_iterations = "\n"

copy = lambda x: bindery.parse(x.xml_encode(), standalone = True)



def read_list (file):
    fh = open(file, mode='r')
    d = {}
    for line in fh:
        d[line.rstrip()]=1
    return d

# parse the XML output and return all query/hit/HSP triplets
def parse(lst):
    # for each query,
    for query in i.BlastOutput.BlastOutput_iterations.Iteration:
        if str(object=query.Iteration_query_ID) not in lst:
            continue
        hits = query.Iteration_hits.Hit
        if (hits is None):
            continue

        # create a copy of this query alone (without children hits)
        query_ = copy(query).Iteration
        query_.Iteration_hits = "\n"

        # for each hit,
        for hit in hits:
            hsps = hit.Hit_hsps.Hsp
            if (hsps is None):
                continue

            # create a copy of this hit alone (without children HSPs)
            hit_ = copy(hit).Hit
            hit_.Hit_hsps = "\n"

            # for each HSP,
            for hsp in hsps:
                yield query_, hit_, hsp

# cast XML tag values into proper type
def cast (value):
    try:
        value = float(value) # test if this is a float
        if (value == int(value)): # test if this is an integer
            value = int(value)
    except:
        pass

    return str(value)

# process the triplets
def process(out,lst):
    has_content, has_query, has_hit = {}, {}, {}
    triplets_per_fn = {}
    fn = os.path.abspath(out)

    for (query, hit, hsp) in parse(lst):
        query_key = str(query.Iteration_iter_num)
        hit_key = str(hit.Hit_num)

        if (fn not in triplets_per_fn):
            triplets_per_fn[fn] = 0

        triplets_per_fn[fn] += 1

        # has anything been added in this output file yet?
        if (fn not in has_content):
            # if not, we create a new document as a copy
            # of the 'envelop' document created earlier
            o = copy(envelop)
            has_content[fn] = o
            has_query[fn] = {}
            has_hit[fn] = {}
        else:
            o = has_content[fn]

        # has this query been added in this output file yet?
        if (query_key not in has_query[fn]):
            # if not, we add to the document the
            # description of the query sequence
            query_ = copy(query).Iteration
            has_content[fn].BlastOutput.BlastOutput_iterations.xml_append(query_)
            has_query[fn][query_key] = query_
            has_hit[fn][query_key] = {}

        # has this hit been added in this output file yet?
        if (hit_key not in has_hit[fn][query_key]):
            # if not, we add to the document the
            # description of the hit sequence
            hit_ = copy(hit).Hit
            has_query[fn][query_key].Iteration_hits.xml_append(hit_)
            has_hit[fn][query_key][hit_key] = hit_

        # we add this HSP to the output file
        hsp_ = copy(hsp).Hsp
        has_hit[fn][query_key][hit_key].Hit_hsps.xml_append(hsp_)
        # yield query_key, hit_key, hsp_key, triplets_per_fn

    # write the output files
    for fn in sorted(has_content):
        has_content[fn].xml_normalize()
        has_content[fn].xml_write(
            writer = amara.writers.lookup("xml-indent"),
            stream = open(fn, 'w')
        )

fof_fh = open(p.fof, mode='r')
for line in fof_fh:
    array = line.split("\t")
    print('treating ' + array[0] + '...')
    qry_lst = read_list(array[0])
    process(array[1].rstrip(),qry_lst)
