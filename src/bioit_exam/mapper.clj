(ns bioit-exam.mapper
  (:require [bioit-exam.spec :as domain]
            [clojure.algo.generic.functor :as f :only fmap]
            [clojure.core.memoize :as memo :only memo]
            [clojure.core.reducers :as r]
            [clojure.spec.alpha :as s]))

(defn compute-kmers
  "Returns all unique kmers with length k mapped to the
  starting indices at which they occur in refseq."
  [refseq k]
  (->> refseq
       (partition k 1)
       (map-indexed (fn [i kmer] {:i i :kmer kmer}))
       (group-by :kmer)
       (f/fmap (partial map :i))))

(def kmer->indices (memo/memo compute-kmers))

(defn get-kmer+indices
  "Returns the indices at which kmer occurs in refseq"
  [refseq kmer]
  (let [indices (kmer->indices refseq (count kmer))]
    (get indices kmer)))

(defn pad-with [fill-in coll]
  (concat coll (repeat fill-in)))

(defn diff-count [coll other]
  (->> (pad-with nil other)
       (map = coll)
       (filter false?)
       (count)))

(defn tolerated?
  "Returns true if readseq has max-miss mismatches at most when
  comparing to refseq at index i."
  [refseq readseq max-miss i]
  (->> (+ i (count readseq))
       (min (dec (count refseq)))
       (subvec refseq i)
       (diff-count readseq)
       (>= max-miss)))

(defn match-indices
  "Find all indices in the given refseq at which a kmer of
  length k matches a kmer in readseq with tolerance of max-miss
  mismatches.
  Be aware that the seed of length k has to match exactly on
  any kmer in readseq, otherwise the whole read is ignored."
  [k max-miss refseq readseq]
  (->> (take k readseq)
       (get-kmer+indices refseq)
       (r/filter #(tolerated? refseq readseq max-miss %))))

(defn- group-by-match-index [indices]
  (fn [mapping readseq]
    (r/reduce #(update %1 %2 conj readseq)
              mapping
              (indices readseq))))

(defn concat-merge
  ([] {})
  ([coll other] (merge-with concat coll other)))

(defn map-reads
  "Maps all reads to refseq and returns a map of indices to the
  reads with a match to refseq at that index."
  [k max-miss {refseq :seq refname :name} reads]
  (let [refvec (vec refseq)
        readseqs (map :seq reads)
        indices #(match-indices k max-miss refvec %)]
    {:reference {:name refname, :seq refvec}
     :mapped-reads (r/fold concat-merge
                           (group-by-match-index indices)
                           readseqs)}))

;; Following is not needed
;; just experimenting a bit

(s/fdef map-reads
  :args (s/cat :k pos-int?
               :max-miss nat-int?
               :refseq ::domain/reference
               :reads (s/coll-of :gen/base))
  :ret ::domain/mapping)

(s/fdef match-indices
  :args (s/cat :k pos-int?
               :max-miss nat-int?
               :refseq ::domain/reference
               :readseq ::domain/read)
  :ret (s/coll-of nat-int?)
  :fn (s/and (s/coll-of nat-int?)
             (fn [x] (every? #(< % (count (-> x :args :refseq))) (:ret x)))))

;; PLAYGROUND

(comment
; (stest/instrument `map-reads)
; (stest/instrument `match-indices)

;; NOTE: tolerance relates to the whole read, not to the kmer!

  (def reference
    {:name "REF001"
     :seq [\X \Y \Z \A \B \C \D \F \F \F \G \A \B \C \D \E \F \D \E \F \F]})
;                  ^  ^  ^  ^  x  x  x              ^  ^  ^  x  ^
;                                          ^  ^  ^  ^  ^  x  ^

  (map-reads 3 1 reference ["ABCDEED" "DEFFE" "GDFE"])

  (match-indices 4 1 "ABCDEFG" "ABCDDEFG"))
