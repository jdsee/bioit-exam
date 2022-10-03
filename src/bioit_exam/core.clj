(ns bioit-exam.core
  (:require [bioit-exam.mapper :refer :all]
            [bioit-exam.polisher :refer :all]
            [bioit-exam.reader :refer :all]
            [bioit-exam.writer :refer :all]))

(defn -main
  [ref-file reads-file k max-miss min-freq]
  (let [reference (first (read-fasta! ref-file))
        reads (read-fasta! reads-file)]
    (->> reads
         ; (polish-all k min-freq)
         (map-reads k max-miss reference)
         (as-pileup)
         (write-pileup! "target/out/result.pileup"))))

(let [k 15, max-miss 2, min-freq 10]
  (-main "data/fluA.fasta"
         "data/fluA_reads.fasta"
         k
         max-miss
         min-freq))
