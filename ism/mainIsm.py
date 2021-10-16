
# MAIN FUNCTION TO CALL THE ISM MODULE

from ism.src.ism import ism

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = '/home/luss/EODP/EODP21/auxiliary'
indir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/' # small scene
outdir = '/home/luss/my_shared_folder/test_ISM/'

# Initialise the ISM
myIsm = ism(auxdir, indir, outdir)
myIsm.processModule()
