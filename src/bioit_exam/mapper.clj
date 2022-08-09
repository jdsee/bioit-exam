(ns bioit-exam.mapper
  (:require [clojure.algo.generic.functor :as f :only fmap]))

(defn read-fasta [file-name])
  ; (with-open [r (BufferedReader. (FileReader.  file-name))]
  ; (.read )))

; (defn map-from-file [ref-file sample-file]
;   (let [reference (read-fasta ref-file)
;         reads (read-fasta sample-file)])
  ; (read-fasta ref-file)
  ; (read-fasta sample-file)

(defn compute-kmers [ref' kmer-size]
  (println "working hard, all day long...")
  (->> ref'
       (partition kmer-size 1)
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

(defn count-mismatches [ref' read' offset]
  (->> ref'
       (drop offset)
       (take (count read'))
       (pad-with nil)
       (map = read')
       (filter false?)
       (count)))

(defn map-reads [kmer-size max-mismatches ref' reads]
  (let [mismatches-tolerated? #(>= max-mismatches (count-mismatches ref' %1 %2))
        matching-positions #(->> (take kmer-size %)
                                 (get-kmer-positions ref')
                                 (filter (partial mismatches-tolerated? %)))]
    (zipmap reads (map matching-positions reads))))

(map-reads 2 2 "cagtcag" ["ag" "gtc" "cag" "aaa"])

