#!/usr/bin/python3

from cobraa.load import *
from cobraa.io_operations import *
from cobraa.coincidence import *
from cobraa.sensitivity import *
from cobraa.globals import *
from cobraa.extras import *

######################## Start of main function ###########################

if __name__ == '__main__':

    print(docstring)
    if arguments['-v']:
    	print(arguments)
#     print defaultValues
    print('')

    if arguments['-m']:
        generateMacros()

    if arguments['-j']:
        generateJobs()

    if arguments['-M']:
        mergeRootFiles()

    if arguments['--coincidences']:
        coincidenceMap()

    if arguments['--sensitivity']:
        calculateSensitivity()

    if arguments['--triggers']:
        triggers()

    print('''All done.

You can submit multiple jobs with a single command line 'for i in `ls job/job*.sh`; do bsub $i; done')

You can use the stand-along version of bonsai (bonsai2) by copying a like.bin to the base folder and doing the following steps:
>bsub
>source /p/gpfs1/adg/wmutils/env.sh
>for i in `ls root_files_*/*/*`; do bonsai2 $i bonsai_$i 3000 800 -500 1000 1;done
>Ctrl-D
which will work on lassen. (Option should be 'do bonsai2 $i bonsai_$i 3000 800 -500 1000 0' for wbls.)
''')

######################## End of main function  ###########################
