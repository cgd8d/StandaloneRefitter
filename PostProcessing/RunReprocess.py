"""
Generate lots of EXOs based on Reprocess.exo.
Then, submit batch jobs to execute all of them.
"""

from __future__ import with_statement
import os, glob, string, subprocess

basepath = '/nfs/slac/g/exo_data4/users/cgd8d/rhel5-64'

os.environ['INPUT_LOCATION'] = os.path.join(basepath, 'DenoisedFiles')
os.environ['OUTPUT_LOCATION'] = os.path.join(basepath, 'ReprocessedFiles')
with open(os.path.join(basepath, 'Refitter', 'PostProcessing', 'Reprocess.exo'), 'r') as InEXO:
    EXOTemplate = string.Template(InEXO.read())

for filename in glob.glob(os.path.join(basepath, 'DenoisedFiles', '*')):
    basefile = os.path.basename(filename)
    os.environ['FILENAME'] = basefile
    with open(os.path.join(basepath, 'EXOs', basefile[:-5]+'.exo'), 'w') as OutEXO:
        OutEXO.write(EXOTemplate.substitute(os.environ))

for filename in glob.glob(os.path.join(basepath, 'EXOs', '*')):
    subprocess.call(['bsub', '-o', filename[:-3]+'log', '-q', 'long', '-R', 'rhel50', 'EXOAnalysis', filename])

