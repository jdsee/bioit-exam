(ns bioit-exam.writer
  (:require
   [clojure.java.io :as io]
   [clojure.pprint :refer :all]
   [clojure.spec.alpha :as s]
   [clojure.spec.test.alpha :as stest]
   [clojure.string :as str]))

(def mapping
  {:ref {:seq [\T \C \A \G \T \G \G \T \C \C \G \T] , :name "REF001"}
   :mapped-reads {5 (list "GGTATG")
              2 (list "AGTAAT" "AGTAATCCGTAG")}})

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
     (if-some [x tail]
       (conj remainder x)
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

(make-pileup mapping)

(s/fdef make-pileup
  :args :core/mapping
  :ret boolean?)

(defn write-pileup [filename pileup]
  (with-open [w (io/writer filename)]
    (for [p pileup
          :let [line (str/join \tab p)]]
      (.write w line))))

(str/join \tab "HALLO")

(stest/instrument `make-pileup)
