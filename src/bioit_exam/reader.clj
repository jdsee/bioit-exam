(ns bioit-exam.reader
  (:require [clojure.java.io :as io]
            [clojure.string :as str]))

(defn header? [line] (= \> (first line)))

(defn last-index [coll] (-> coll count dec))

(defn extr-name
  [header]
  (second (re-find #">\s*\b(\S+)\b.*" header)))

(defn parse-fasta-line
  [reads line]
  (let [trimmed (str/trim line)]
    (if (header? trimmed)
      (conj reads {:header line
                   :name (extr-name trimmed)
                   :seq ""})
      (update-in reads
                 [(last-index reads) :seq]
                 #(str % trimmed)))))

(defn parse-fasta
  [lines]
  (->> (remove empty? lines)
       (reduce parse-fasta-line [])))

(defn read-fasta! [file-name]
  (with-open [r (io/reader file-name)]
    (parse-fasta (line-seq r))))

(read-fasta! "data/fluA_reads.fasta")
