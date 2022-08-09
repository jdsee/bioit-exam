from functools import reduce
from mapperutil import *
import operator
import re

class Sequence:
    def __init__(self, lines):
        maybe_name = re.search(r">\s*\b(\S+)\b.*", lines[0])
        self.name = maybe_name.group(1) if maybe_name else "unknown"
        self.bases = ''.join(lines[1:])

    def __str__(self):
        dots = '...' if len(self.bases) >= 20 else ''
        return f'{self.name}: {self.bases[0:20]}{dots}'

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.bases)


class Read(Sequence):
    def get_seed(self, seedlength):
        return self.bases[:seedlength]

    def replace_kmers(self, replacements):
        for kmer, rep in replacements:
            for err in re.finditer(kmer, self.bases):
                self.bases = replace(self.bases, rep, err.start(), err.end())


class Reference(Sequence):
    def __init__(self, lines):
        self.kmersize = None
        self.kmers = {}
        super().__init__(lines)

    def calculate_kmers(self, kmersize):
        self.kmersize = kmersize
        for i in range(len(self.bases) - kmersize + 1):
            kmer = self.bases[i:i+kmersize]
            self.kmers[kmer] = self.kmers.get(kmer, []) + [i]

    def get_kmer_positions(self, kmer):
        seedlen = len(kmer)
        if (self.kmersize != seedlen):
            self.calculate_kmers(seedlen)
        return self.kmers.get(kmer, [])

    def count_mismatches(self, read: Read, position: int):
        snippet  = self.bases[ position : position+len(read.bases) ]
        overflow = len(read.bases) - len(snippet)
        return sum(map(operator.ne, read.bases, snippet)) + overflow


class Mapping:
    def __init__(self, reference):
        self.ref = reference
        self.reads = {}

    def add_read(self, read, position):
        self.reads[position] = self.reads.get(position, []) + [read]

    def get_reads_at_position(self, position):
        return self.reads.get(position, [])

    def get_pileup(self):
        res, rem = [], []
        for i, base in enumerate(self.ref.bases):
            rem += [r.bases for r in self.get_reads_at_position(i)]
            matches, tmp = "", []
            for read in rem:
                matches += "." if read[0:1] == base else read[0:1]
                tmp.append(read[1:])
            rem = tmp
            res.append([i + 1, base, len(matches), matches])
        return res

    def __str__(self):
        res = ["Mapping to " + self.ref.name]
        for pos in self.reads:
            res += ["  " + str(len(self.reads[pos])) + " reads mapping at " + str(pos)]
        return "\n".join(res)

    def __repr__(self) -> str:
        return str(self)


class ReadPolisher:
    def __init__(self, kmerlen):
        self.kmerlen = kmerlen
        self.spectrum = {}

    def add_read(self, read):
        last_window_start = len(read.bases) - self.kmerlen + 1
        for i in range(last_window_start):
            kmer = read.bases[i : i + self.kmerlen]
            count, reads = self.spectrum.get(kmer, (0, set()))
            reads.add(read)
            self.spectrum[kmer] = count + 1, reads

    def get_replacements(self, minfreq):
        candidates = self._get_candidates(minfreq)
        reps = self._filter_replacements(minfreq, candidates)
        return groupby_key(reps)

    def _filter_replacements(self, minfreq, candidates):
        by_count  = lambda x : self._get_count(x[0])
        acc_reads = lambda acc, x : acc | x[1]
        return (
            (read, (kmer, rep))
            for kmer, rep_reads in groupby_key(candidates).items()
            if  (best := max(rep_reads, key=by_count))
            and self._get_count(rep := best[0]) >= minfreq
            for read in reduce(acc_reads, rep_reads, set())
        )

    def _get_candidates(self, minfreq):
        return (
            (kmer, (replace(kmer, base, i), reads))
            for kmer, (count, reads) in self.spectrum.items()
            if  count < minfreq
            for i in range(len(kmer))
            for base in "AGTC"
            if  base != kmer[i]
        )

    def _get_count(self, kmer):
        count, _ = self.spectrum.get(kmer, (0, None))
        return count


def read_fasta(readfile, classname):
    reads = []
    with open(readfile) as file:
        for line in file:
            stripped = line.strip()
            if stripped and stripped[0] == '>':
                reads.append([stripped])
            else:
                reads[-1].append(stripped)
        return list(map(globals()[classname], reads))

def map_reads(reads, reference, kmersize, max_mismatches):
    mapping = Mapping(reference)
    for read in reads:
        seed = read.get_seed(kmersize)
        positions = reference.get_kmer_positions(seed)
        for pos in positions:
           mismatches = reference.count_mismatches(read, pos)
           if mismatches < max_mismatches:
               mapping.add_read(read, pos)
    return mapping


class MappingWriter:
    def __init__(self, mapping):
        self.mapping = mapping

    def write_sam(self, filename):
        ref = self.mapping.ref
        sq = f"@SQ\tSN:{ref.name}\tLN:{len(ref.bases)}"
        mappings = "\n".join(\
            self._sam_line(
                qname=read.name,
                rname=ref.name,
                pos=pos + 1,
                cigar=f"{len(read)}M",
                seq=read.bases)
            for pos, reads in self.mapping.reads.items()
            for read in reads)
        with open(filename, 'w') as file:
            file.write(f"{sq}\n{mappings}")
            file.write(f"\n{mappings}")

    def write_pileup(self, filename):
        with open(filename, 'w') as file:
            name = self.mapping.ref.name
            pileups = self.mapping.get_pileup()
            for pileup in pileups:
                file.write('{}\t{}\t{}\t{}\t{}\n'.format(name, *pileup))

    def _sam_line(self, qname, rname, pos, cigar, seq, flag = 0, mapq = 255, rnext = "*", pnext = 0, tlen = 0, qual = "*"):
        return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"


def main():
    reads = read_fasta("data/patient1.fasta", Read.__name__)
    reference = read_fasta("data/rpoB.fasta", Reference.__name__)[0]
    mapping = map_reads(reads, reference, 158, 2)
    writer = MappingWriter(mapping)
    writer.write_sam("data/mapping_p1_uncorrected.sam")
    writer.write_pileup("data/mapping_p1_uncorrected.pileup")
    polisher = ReadPolisher(15)
    for read in reads:
        polisher.add_read(read)
    replacements = polisher.get_replacements(10)
    for read, reps in replacements.items():
        read.replace_kmers(reps)
    mapping = map_reads(reads, reference, 15, 2)
    writer = MappingWriter(mapping)
    writer.write_sam("data/mapping_p1_corrected.sam")
    writer.write_pileup("data/mapping_p1_corrected.pileup")


if __name__ == "__main__":
    main()

