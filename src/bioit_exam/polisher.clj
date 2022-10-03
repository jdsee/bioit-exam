(ns bioit-exam.polisher
  (:require [clojure.algo.generic.functor :as f]))

(defn slice-kmers [k read']
  (->> read'
       (partition k 1)
       (map (fn [kmer] {:kmer kmer :read read'}))))

(defn build-spectrum [k reads]
  (->> reads
       (mapcat (partial slice-kmers k))
       (group-by :kmer)
       (f/fmap #(let [reads (map :read %)]
                  {:reads (distinct reads)
                   :freq (count reads)}))))

(defn polish [kmer replacement pos]
  (let [head (take pos kmer)
        tail (drop (inc pos) kmer)]
    (concat head replacement tail)))

(defn candidates [min-freq spectrum]
  (for [[kmer {freq :freq, reads :reads}] spectrum
        i (range (count kmer))
        base [\A \G \T \C]
        :when (and (< freq min-freq) (not= base (get kmer i)))]
    {:kmer kmer
     :replacement (polish kmer base i)
     :reads reads}))
;; TODO: remove where (= :kmer :replacement)

(defn best-candidate [replace-freqs candidates]
  (apply max-key replace-freqs candidates))

(defn choose-candidates [min-freq spectrum]
  (let [replace-freqs #(get-in spectrum [(:replacement %) :freq] 0)]
    (->> (candidates min-freq spectrum)
         (group-by :kmer)
         (vals)
         (map (partial best-candidate replace-freqs))
         (filter #(>= (replace-freqs %) min-freq)))))

(defn polish-all
  [k min-freq reads]
  (->> (build-spectrum k reads)
       (choose-candidates min-freq)))

