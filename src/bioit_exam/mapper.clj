(ns bioit-exam.mapper
  (:require [clojure.algo.generic.functor :as f :only fmap]
            [clojure.core.memoize :as memo :only memo]
            [clojure.spec.alpha :as s]
            [clojure.spec.test.alpha :as stest]))

(defn compute-kmers
  "Returns all unique kmers with length k mapped to the
  starting postions at which they occur in refseq."
  [refseq k]
  (->> refseq
       (partition k 1)
       (map-indexed (fn [i kmer] {:pos i :kmer kmer}))
       (group-by :kmer)
       (f/fmap (partial map :pos))))

(def kmer-positions (memo/memo compute-kmers))

(defn get-kmer-positions
  "Returns the positions at which kmer occurs in refseq"
  [refseq kmer]
  (-> (kmer-positions refseq (count kmer))
      (get kmer [])))

(defn pad-with [fill-in coll]
  (concat coll (repeat fill-in)))

(defn diff-count [coll other]
  (->> (pad-with nil other)
       (map = coll)
       (filter false?)
       (count)))

(defn tolerated?
  "Returns true if readseq has max-div mismatches at most when
  comparing to refseq at position i."
  [refseq readseq i max-div]
  (->> (+ i (count readseq))
       (subs refseq i)
       (diff-count readseq)
       (#(<= % max-div))))

(defn matching-positions
  "Find all positions in the given refseq at which kmers of
  length k match kmers in readseq with tolerance of max-div
  mismatches.
  Be aware that the seed of length k has to match exactly on
  any kmer in readseq, otherwise the whole read is ignored."
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

(defn all-caps? [string]
  (every? #(Character/isUpperCase %) string))

(s/def :gen/base (s/and string? all-caps?))

(s/fdef map-reads
  :args (s/cat :k pos-int?
               :max-div nat-int?
               :refseq :gen/base
               :reads (s/coll-of :gen/base)
  :ret (s/map-of int? string?)))

(s/fdef matching-positions
  :args (s/cat :k pos-int?
               :max-div nat-int?
               :refseq string?
               :readseq string?)
  :ret (s/coll-of nat-int?)
  :fn (s/and (s/coll-of nat-int?)
             (fn [x] (every? #(< % (count (-> x :args :refseq))) (:ret x)))))
