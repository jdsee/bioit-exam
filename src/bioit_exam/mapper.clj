(ns bioit-exam.mapper
  (:require [clojure.algo.generic.functor :as f :only fmap]
            [clojure.spec.alpha :as s]))

(defn compute-kmers
  "Finds all unique kmers with length k and
  maps them to matching postions in ref-seq."
  [ref-seq k]
  (->> ref-seq
       (partition k 1)
       (map-indexed (fn [i kmer] {:pos i :kmer kmer}))
       (group-by :kmer)
       (f/fmap (partial map :pos))))

(def kmer-positions
  (memoize compute-kmers))

(defn get-kmer-positions [ref-seq seed]
  (-> ref-seq
      (kmer-positions (count seed))
      (get seed [])))

(defn pad-with [fill-in coll]
  (concat coll (repeat fill-in)))

(defn diff-count [coll other]
  (->> (pad-with nil coll)
       (map = (pad-with nil other))
       (filter false?)
       (count)))

(defn tolerated?
  "Return true if read-seq has max-div mismatches at max when
  comparing to ref-seq at position pos."
  [ref-seq read-seq pos max-div]
  (->> (drop pos ref-seq)
       (take (count read-seq))
       (diff-count read-seq)
       #(> % max-div)))

(defn matching-positions
  "Find all positions in the given ref-seq at which kmers of
  length k match kmers in read-seq with tolerance of max-div
  mismatches."
  [k max-div ref-seq read-seq]
  (->> (take k read-seq)
       (get-kmer-positions ref-seq)
       (filter #(tolerated? ref-seq read-seq max-div %))))

(defn group-by-pos [k max-div ref-seq]
  (fn [mapping read-seq]
    (->> (matching-positions k max-div ref-seq read-seq)
         (reduce #(update %1 %2 concat [read-seq]) mapping))))
         ; (map #(update mapping % concat [read-seq])))))

(defn map-reads
  "Maps all reads to ref-seq and returns a map of indices to the
  reads with potential match to ref-seq at that index."
  [k max-div ref-seq reads]
  (reduce (group-by-pos k max-div ref-seq) {} reads))

;; Following is not needed
;; just experimenting a bit

(s/fdef map-reads
  :args (s/cat :k pos-int?
               :max-div nat-int?
               :ref-seq string?
               :reads (s/coll-of string?))
  :ret (s/map-of int? string?))

(s/fdef matching-positions
  :args (s/cat :k pos-int?
               :max-div nat-int?
               :ref-seq string?
               :read-seq string?)
  :ret (s/coll-of nat-int?)
  :fn (s/and (s/coll-of nat-int?)
             (fn [x] (every? #(< % (count (-> x :args :ref-seq))) (:ret x)))))
