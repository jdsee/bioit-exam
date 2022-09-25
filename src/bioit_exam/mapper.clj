(ns bioit-exam.mapper
  (:require [clojure.algo.generic.functor :as f :only fmap]
            [clojure.spec.alpha :as s]))

(defn compute-kmers
  "Returns all unique kmers with length k mapped to the
  starting postions at which they occur in refseq."
  [refseq k]
  (->> refseq
       (partition k 1)
       (map-indexed (fn [i kmer] {:pos i :kmer kmer}))
       (group-by :kmer)
       (f/fmap (partial map :pos))))

(def kmer-positions (memoize compute-kmers))

(defn get-kmer-positions
  "Return the positions at which kmer occurs in refseq"
  [refseq kmer]
  (-> (kmer-positions refseq (count kmer))
      (get kmer [])))

(defn pad-with [fill-in coll]
  (concat coll (repeat fill-in)))

(defn diff-count [coll other]
  (->> (pad-with nil coll)
       (map = (pad-with nil other))
       (filter false?)
       (count)))

(defn tolerated?
  "Return true if readseq has max-div mismatches at most when
  comparing to refseq at position i."
  [refseq readseq i max-div]
  (->> (subvec refseq i (count readseq))
       (diff-count readseq)
       #(> % max-div)))

(defn matching-positions
  "Find all positions in the given refseq at which kmers of
  length k match kmers in readseq with tolerance of max-div
  mismatches."
  [k max-div refseq readseq]
  (->> (take k readseq)
       (get-kmer-positions refseq)
       (filter #(tolerated? refseq readseq max-div %))))

(defn group-by-pos [k max-div refseq]
  (fn [mapping readseq]
    (->> (matching-positions k max-div refseq readseq)
         (reduce #(update %1 %2 concat [readseq]) mapping))))

(defn map-reads
  "Maps all reads to refseq and returns a map of indices to the
  reads with a match to refseq at that index."
  [k max-div refseq reads]
  (reduce (group-by-pos k max-div refseq) {} reads))

;; Following is not needed
;; just experimenting a bit

(s/fdef map-reads
  :args (s/cat :k pos-int?
               :max-div nat-int?
               :refseq string?
               :reads (s/coll-of string?))
  :ret (s/map-of int? string?))

(s/fdef matching-positions
  :args (s/cat :k pos-int?
               :max-div nat-int?
               :refseq string?
               :readseq string?)
  :ret (s/coll-of nat-int?)
  :fn (s/and (s/coll-of nat-int?)
             (fn [x] (every? #(< % (count (-> x :args :refseq))) (:ret x)))))
