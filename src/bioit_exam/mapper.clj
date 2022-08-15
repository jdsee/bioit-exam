(ns bioit-exam.mapper
  (:require [bioit-exam.polisher :refer :all]
            [clojure.algo.generic.functor :as f :only fmap]))

(defn compute-kmers [ref' kmerlen]
  (->> ref'
       (partition kmerlen 1)
       (map-indexed (fn [i kmer] {:pos i :kmer kmer}))
       (group-by :kmer)
       (f/fmap (partial map :pos))))

(def kmer-positions
  (memoize compute-kmers))

(defn get-kmer-positions [ref' kmer]
  (-> ref'
      (kmer-positions (count kmer))
      (get kmer [])))

(defn pad-with [fill-in xs]  (concat xs (repeat fill-in)))

(defn count-mismatches [ref' read' pos]
  (->> ref'
       (drop pos)
       (take (count read'))
       (pad-with nil)
       (map = read')
       (filter false?)
       (count)))

(defn matching-positions [kmerlen max-mismatches ref' read']
  (->> read'
       (take kmerlen)
       (get-kmer-positions ref')
       (filter #(>= max-mismatches (count-mismatches ref' read' %)))))

(defn map-reads [kmerlen max-mismatches ref' reads]
  (->> reads
       (map (partial matching-positions kmerlen max-mismatches ref'))
       (zipmap reads)))

