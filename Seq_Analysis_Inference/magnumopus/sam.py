import re
from collections import defaultdict
from copy import deepcopy


                            ##############################
                            ##### Read + SAM classes #####
                            ##############################

class Read:
	def __init__(self, sam_line: str):
		field = sam_line.strip().split("\t") 

        #store each field of the SAM entry in named attributes 
		self.qname = field[0] 
		self.flag = int(field[1])
		self.rname = field[2]
		self.pos = int(field[3])
		self.mapq = int(field[4])
		self.cigar = field[5]
		self.rnext = field[6]
		self.pnext = int(field[7])
		self.tlen = int(field[8])
		self.seq = field[9]
		self.qual = field[10]

        #mapping properties based on flags
		self.is_mapped = not bool(self.flag & 4) #4 bit not in flag
		self.is_forward = not bool(self.flag & 16)
		self.is_reverse = bool(self.flag & 16)
		self.is_primary = not (bool(self.flag & 256) or bool(self.flag & 2048)) #not secondary and not supplemental


		#add data for mapped reads only
		self.cigar_bits: tuple[tuple[int, str]] = None
		self.mapped_len: int = None
		
		if self.is_mapped:
			#parse cigar string
			self.cigar_bits = tuple([(int(n), cig) for n, cig in re.findall(r"(\d+)([A-Z])", self.cigar)])
			#find mapped length
			self.mapped_len = sum([n for n, cig in self.cigar_bits if cig in {"M", "D"}])


	def read_idx_at_pos(self, pos: int) -> list[None|int]:
		"""
            takes a 1-base position in the reference sequence and returns the corresponding read index
        """
		
		#if not mapped, exit
		if not self.is_mapped:
			return []
		
		# adjust by read start
		pos -= self.pos
		if pos < 0: # If read mapped to the right of requested location
			return []

		# Check if the requested position is right of our read
		if pos >= self.mapped_len:
			return []
		
		cigar_bits = re.findall(r"(\d+)([A-Z])", self.cigar)
		# seq position based on cigar
		seq = 0
		mapped_count = 0
		for n, (size, cig_type) in enumerate(cigar_bits):
			size = int(size)
			if cig_type in {"S", "H", "I"}: #soft clippings, hard clippings, insertions (only read consumed)
				seq += size
				continue
			if mapped_count + size >= pos+1:
				if cig_type == "M": #match/mismatch
					# seq remaining
					seq += (pos-mapped_count)
					break
				if cig_type == "D": #deletion
					return []
			else:
				if cig_type == "M":
					mapped_count += size
					seq += size
				if cig_type == "D":
					mapped_count += size

		# Check if next bases are insertion
		if mapped_count + size == pos+1:
			if n+1 != len(cigar_bits):
				size, cig_type = cigar_bits[n+1]
				size = int(size)
				if cig_type == "I":
					return [i for i in range(seq, seq+size+1)]
		
		return [seq]


	def mapped_seq(self) -> str:
		"""
        returns the mapped portion of the read that corresponds to the reference sequence
        """
		#if not mapped, exit
		if not self.is_mapped:
			return ""

		idx = 0 # track where we are in read seq 
		bases = [] # list to build up over time
		
		for len, cig in self.cigar_bits:
			if cig == "S":
				idx += len
			elif cig == "D":
				bases += ["-"]*len
			elif cig in {"M", "I"}:
				bases += [self.seq[i] for i in range(idx, idx+len)]
				idx += len

		return "".join(bases)


	def base_at_pos(self, pos: int) -> str:
		"""
            takes a 1-base position in the reference sequence and returns the base mapped to specific position in the reference (same strand)
        """
		idx = self.read_idx_at_pos(pos)
		return "".join([self.seq[i] for i in idx])

	
	def qual_at_pos(self, pos: int) -> str:
		"""
        takes a 1-base position in the reference sequence and returns the quality score of a base mapped to a specific position in the reference
        """
		idx = self.read_idx_at_pos(pos)
		return "".join([self.qual[i] for i in idx])


class SAM:
	def __init__(self):
		self._reads: dict[str, dict[int, list[Read]]] = defaultdict(lambda: defaultdict(list))
	

	@property
	def reads(self):
		return deepcopy(self._reads)
	
	
	@property
	def best_ref(self):
		"""
		returns name of best reference
		"""
		most_reads = 0
		best_ref = ""
		for ref_name, pos in self.reads.items():
			if len(pos) > most_reads:
				most_reads = len(pos)
				best_ref = ref_name
		return best_ref


	def add_read(self, read: Read):
		"""
		add read to reference position
		"""
		if read.is_mapped and read.is_primary:
			for ref_pos in range(read.pos, read.pos + read.mapped_len):
				self._reads[read.rname][ref_pos].append(read)
		else:
			# don't store unmapped reads
			pass


	def reads_at_pos(self, seq_name: str, pos: int):
		"""
        returns all Read instances that map to that position in the reference sequence
        """
		return self._reads[seq_name].get(pos, [])


	def pileup_at_pos(self, seq_name: str, pos: int) -> tuple[list[str], list[str]]:
		"""
        returns the base calls and their corresponding quality scores at a given position in the reference
        """
		s_pileup = []
		q_pileup = []
		#call reads class to get bases and quality scores at input position
		for read in self._reads[seq_name].get(pos, []):
			s_pileup.append(read.base_at_pos(pos))
			q_pileup.append(read.qual_at_pos(pos))
		return (s_pileup, q_pileup)

	
	def consensus_at_pos(self, seq_name: str, pos: int) -> str:
		"""
        returns the majority base call at a position in the reference
        """
		#quit if unmapped
		if pos not in self._reads[seq_name]:
			return ""
		
		#get bases and qualities
		bases, quals = self.pileup_at_pos(seq_name, pos)

		#count bases
		base_count = defaultdict(int)
		for base, qual in zip(bases, quals):
			# Apply qual cutoff here if needed
			base_count[base] += 1
		
		#sort bases by count
		sorted_base_count = sorted([(b, c) for b, c in base_count.items()], key=lambda x: x[1], reverse=True)
		if len(sorted_base_count) > 1:
			if sorted_base_count[0][1] <= sum([i[1] for i in sorted_base_count])/2:
				# 50% or more bases suggest this is not the base
				top_base = "N"
				return top_base

		top_base = sorted_base_count[0][0]
		return top_base


	def consensus(self, seq_name: str, start: int=None, end: int=None, fasta=False):
		"""
        returns the majority base call for all positions in the specified reference 
        """
		if fasta:
			consensus = f">{seq_name}_consensus\n"
		else:
			consensus = ""

		if seq_name not in self.reads:
			if fasta:
				consensus += "\n"
			return consensus

		if start is None or start < 1:
			start = 1
		if end is None or end > max(self._reads[seq_name].keys()):
			end = max(self._reads[seq_name].keys())
		if fasta:
			consensus = f">{seq_name}_consensus\n"
		consensus = "".join([consensus] + [self.consensus_at_pos(seq_name, i) for i in range(start, end + 1)])
		if fasta:
			consensus += "\n" 
		return consensus
	

	def best_consensus(self, fasta=False):
		"""
        returns the majority base call for all positions in the reference with the best mapping (i.e., the more positions reads map to, the better the mapping)
        """
		best = self.best_ref
		return self.consensus(best, fasta=fasta)


	@classmethod
	def from_sam(cls, sam_file: str):
		"""
        creates an instance of SAM from a SAM file (and Read instance for each entry in SAM file)
        """
		sam = cls()
		with open(sam_file) as fin:
			for line in fin:
				if line[0] != "@":
					# process SAM line
					read = Read(line)
					sam.add_read(read)
					continue
		return sam




