(ns bioit-exam.writer
  (:require
   [clojure.java.io :as io]
   [clojure.string :as str]))

(defn- make-line
  [{refname :name refseq :seq} i mapped-bases]
  {:refname refname
   :pos (inc i)
   :refbase (get refseq i)
   :matches (count mapped-bases)
   :mapped-bases mapped-bases})

(defn chop-matching-bases [refbase coll]
  (reduce
   (fn [[bases stash] [base & remain]]
     [(str bases (str/replace base refbase \.))
      (if (some? remain) (conj stash remain) stash)])
   ["" []]
   coll))

(defn- pileup-bases
  [{:keys [mapped-reads] {refseq :seq} :reference} stash i]
  (->> (mapped-reads i)
       (concat stash)
       (chop-matching-bases (get refseq i))))

(defn make-pileup
  [{{refseq :seq :as reference} :reference :as mapping}]
  (loop [stash [], accum [], i 0]
    (let [[bases stash] (pileup-bases mapping stash i)
          next-pos (inc i)
          accum (conj accum (make-line reference i bases))]
      (if (> (count refseq) next-pos)
        (recur stash accum next-pos)
        accum))))

(defn pileup-lines
  [pileup]
  (map #(str/join \tab (vals %)) pileup))

(defn as-pileup
  [mapping]
  (-> mapping make-pileup pileup-lines))

(defn write-pileup! [filename pileup-lines]
  (with-open [w (io/writer filename)]
    (doseq [line pileup-lines]
      (.write w line)
      (.newLine w))))

; (def mapping
;   {:reference {:seq [\T \C \A \G \T \G \G \T \C \C \G \T] , :name "REF001"}
;    :mapped-reads {5 (list "GGTATG")
;                   2 (list "AGTAAT" "AGTAATCCGTAG")}})

; (->> (make-pileup mapping)
;      (pileup-lines)
;      (write-pileup! "test_pileup"))
