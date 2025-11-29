import tempfile
import subprocess

from magnumopus.sam import SAM


                            #############################
                            ##### mapping functions #####
                            #############################

def minimap(ref: str, read1: str, read2: str) -> str:
    """
    maps two reads to a reference and returns the best consensus sequence 
    """

    minimaps_sp = ["minimap2", "-ax", "sr", "-B", "1","--score-N", "0", "-k", "10", ref, read1, read2]
    result = subprocess.run(minimaps_sp, capture_output=True, text=True)
    sam_output = result.stdout

    # write SAM to temp file
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".sam") as temp:
        temp.write(sam_output)
        sam_path = temp.name

    #create sam instance 
    sam1 = SAM.from_sam(sam_path) 

    consensus = sam1.best_consensus()

    #return sam1
    return consensus
