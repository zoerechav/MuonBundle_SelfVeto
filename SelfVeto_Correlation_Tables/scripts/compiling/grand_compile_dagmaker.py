import os
import sys
import glob
import numpy as np

dataset = '22803'
subdir = '0000000-0004999'

def get_Grid_job(infile,outfile):
    
    job_name = 'grand_compile'
    lines =[
        'JOB ' + job_name + ' /home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/grand_compile.sub',
        'VARS ' + job_name + ' infile="' + str(infile) + '" outfile="' + outfile +'" ',
        'Retry ' + job_name + ' 2',
        ]
    return lines


input_path = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/' + dataset + '/*.hdf5'
output_path = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/' + dataset + '/correlation/'+ dataset + '_' + subdir + '_correlated.hdf5'

input_files = sorted(glob.glob(input_path))


dag_lines = []
dag_lines.extend(
get_Grid_job(
    infile=input_files,
    outfile = output_path,

    )
)
outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/dags/' + dataset + '_' + subdir + '_grand_compile.dag'
#outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/dags/' + dataset + '_' + subdir + '_weight_airshowered_corsika.dag'
with open(outfile_name, 'w') as f:
    f.write('\n'.join(dag_lines))
