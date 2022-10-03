(ns bioit-exam.mapper
  (:require [bioit-exam.spec :as domain]
            [clojure.algo.generic.functor :as f :only fmap]
            [clojure.core.memoize :as memo :only memo]
            [clojure.core.reducers :as r]
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
  (let [positions (kmer-positions refseq (count kmer))]
    (get positions kmer)))

(defn pad-with [fill-in coll]
  (concat coll (repeat fill-in)))

(defn diff-count [coll other]
  (->> (pad-with nil other)
       (map = coll)
       (filter false?)
       (count)))

(defn tolerated?
  "Returns true if readseq has max-miss mismatches at most when
  comparing to refseq at position i."
  [refseq readseq max-miss i]
  (->> (+ i (count readseq))
       (min (dec (count refseq)))
       (subvec refseq i)
       (diff-count readseq)
       (>= max-miss)))

(defn matching-positions
  "Find all positions in the given refseq at which kmers of
  length k match kmers in readseq with tolerance of max-miss
  mismatches.
  Be aware that the seed of length k has to match exactly on
  any kmer in readseq, otherwise the whole read is ignored."
  [k max-miss refseq readseq]
  (->> (take k readseq)
       (get-kmer-positions refseq)
       (r/filter #(tolerated? refseq readseq max-miss %))))

(defn- group-by-pos [positions]
  (fn [mapping readseq]
    (r/reduce #(update %1 %2 conj readseq)
              mapping
              (positions readseq))))

(defn concat-merge
  ([] {})
  ([coll other] (merge-with concat coll other)))

(defn map-reads
  "Maps all reads to refseq and returns a map of indices to the
  reads with a match to refseq at that index."
  [k max-miss {refseq :seq refname :name} reads]
  (let [refvec (vec refseq)
        readseqs (map :seq reads)
        positions #(matching-positions k max-miss refvec %)]
    {:reference {:name refname, :seq refvec}
     :mapped-reads (r/fold concat-merge
                           (group-by-pos positions)
                           readseqs)}))

;; Following is not needed
;; just experimenting a bit

(s/fdef map-reads
  :args (s/cat :k pos-int?
               :max-miss nat-int?
               :refseq :domain/reference
               :reads (s/coll-of :gen/base))
  :ret :domain/mapping)

(s/fdef matching-positions
  :args (s/cat :k pos-int?
               :max-miss nat-int?
               :refseq :domain/reference
               :readseq :domain/read)
  :ret (s/coll-of nat-int?)
  :fn (s/and (s/coll-of nat-int?)
             (fn [x] (every? #(< % (count (-> x :args :refseq))) (:ret x)))))

;; PLAYGROUND

(stest/instrument `map-reads)
; (stest/instrument `matching-positions)

;; NOTE: tolerance relates to the whole read, not to the kmer!

(def reference
  {:name "REF001"
   :seq [\X \Y \Z \A \B \C \D \F \F \F \G \A \B \C \D \E \F \D \E \F \F]})
;                  ^  ^  ^  ^  x  x  x              ^  ^  ^  x  ^
;                                          ^  ^  ^  ^  ^  x  ^

(map-reads 3 1 reference ["ABCDEED" "DEFFE" "GDFE"])

(matching-positions 4 1 "ABCDEFG" "ABCDDEFG")
