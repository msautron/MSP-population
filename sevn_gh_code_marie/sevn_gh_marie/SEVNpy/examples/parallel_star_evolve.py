"""
Simple script showing how is easy to evolve a list of stars in parallel using the Python
multiprocess module
"""

import numpy as np
import multiprocess as mp
import time
from sevnpy.sevn import SEVNmanager,Star


def apply_evolve(star):
    star.evolve()
    return star

if __name__=="__main__":

    #Initialise the SEVN session, with default parameter
    SEVNmanager.init()

    #Get minimum and maximum ZAMS mass and Z
    tables_info = SEVNmanager.tables_info()
    min_zams = tables_info["min_zams"]
    max_zams = tables_info["max_zams"]
    min_Z = tables_info["min_z"]
    max_Z = tables_info["max_z"]

    #Number of stars
    Nstars = 10000

    Star_list = [Star(mzams,z) for mzams,z
                 in zip(np.random.uniform(min_zams,max_zams,Nstars),np.random.uniform(min_Z,max_Z,Nstars))]


    #Serial run
    t1=time.perf_counter()
    Stars_evolved=list(map(apply_evolve,Star_list))
    t2=time.perf_counter()
    print(f"Runtime nproc=1: {t2-t1}",flush=True)

    #Parallel run
    nproc=2
    t1=time.perf_counter()
    with mp.Pool(nproc) as pool:
        Stars_evolved_parallel_a=pool.map(apply_evolve,Star_list)
    t2=time.perf_counter()
    print(f"Runtime nproc={nproc}: {t2-t1}",flush=True)

    nproc=4
    t1=time.perf_counter()
    with mp.Pool(nproc) as pool:
        Stars_evolved_parallel_b=pool.map(apply_evolve,Star_list)
    t2=time.perf_counter()
    print(f"Runtime nproc={nproc}: {t2-t1}",flush=True)

    #Close the SEVN session
    SEVNmanager.close()
