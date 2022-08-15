(ns bioit-exam.reader
  (:require [clojure.java.io :as io]
            [clojure.string :as str]))

(defn parse-fasta-line [reads line]
  (let [last-index (-> reads count dec)
        trimmed (str/trim line)]
    (if (= \> (first trimmed))
      (conj reads {:header line :lines []})
      (update-in reads [last-index :lines] #(conj %1 trimmed)))))

(defn parse-fasta [lines]
  (reduce parse-fasta-line [] lines))

(defn read-fasta [file-name]
  (with-open [rdr (io/reader file-name)]
    (parse-fasta (line-seq rdr))))

