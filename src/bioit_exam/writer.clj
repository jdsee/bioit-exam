(ns bioit-exam.writer
  (:require
   [clojure.java.io :as io]
   [clojure.pprint :refer :all]
   [clojure.spec.alpha :as s]
   [clojure.spec.test.alpha :as stest]
   [clojure.string :as str]))

(defn- base-repr [refbase base]
  (if (= base refbase) \. base))

(defn- make-line
  [{refname :name refseq :seq} i mapped-bases]
  {:refname refname
   :pos (inc i)
   :refbase (get refseq i)
   :matches (count mapped-bases)
   :mapped-bases mapped-bases})

(defn- line-partitioner [refbase]
  (fn [[result remainder] [base & tail]]
    [(str result (base-repr refbase base))
     (if (some? tail)
       (conj remainder tail)
       remainder)]))

(defn- pileup-bases
  [{:keys [mapped-reads] {refseq :seq} :ref} stash i]
  (->> (mapped-reads i)
       (concat stash)
       (reduce (line-partitioner (get refseq i)) ["" []])))

(defn make-pileup
  [{{refseq :seq :as reference} :ref :as mapping}]
  (loop [stash "", accum [], i 0]
    (let [[bases stash] (pileup-bases mapping stash i)
          line (make-line reference i bases)
          next-pos (inc i)
          accum (conj accum line)]
      (if (> (count refseq) next-pos)
        (recur stash accum next-pos)
        accum))))

(defn pileup-lines
  [pileup]
  (map #(str/join \tab (vals %)) pileup))

(defn write-pileup [filename pileup-lines]
  (with-open [w (io/writer filename)]
    (doseq [line pileup-lines]
      (.write w line)
      (.newLine w))))
