(ns bioit-exam.polisher
  (:require
   [clojure.algo.generic.functor :as f]
   [clojure.string :as str]))

(defn- slice-kmers
  [k i readseq]
  (->> readseq
       (partition k 1)
       (map #(hash-map % [i]))))

(defn build-spectrum
  [k reads]
  (->> reads
       (map-indexed #(slice-kmers k %1 %2))
       (flatten)
       (apply merge-with concat)
       (f/fmap #(hash-map :indices (distinct %)
                          :freq (count %)))))

(defn polish-kmer
  [kmer subst i]
  (concat (take i kmer)
          [subst]
          (drop (inc i) kmer)))

(defn- candidates
  [min-freq spectrum]
  (for [[kmer {:keys [freq indices]}] spectrum
        i     (range (count kmer))
        base  (disj #{\A \G \T \C} (nth kmer i))
        :let  [alt (polish-kmer kmer base i)
               altfreq (get-in spectrum [alt :freq] 0)]
        :when (and  (< freq min-freq)
                    (>= altfreq min-freq))]
    {:kmer kmer
     :alt alt
     :freq (get-in spectrum [alt :freq] 0)
     :indices indices}))

(defn- find-alternatives
  [min-freq spectrum]
  (->> (candidates min-freq spectrum)
       (group-by :kmer)
       (vals)
       (map #(apply max-key :freq %))
       ; (filter #(get-in spectrum [(:alt %) :freq]))
       ))

(defn fix-reads
  [reads {:keys [kmer alt indices]}]
  (let [kmerstr (str/join kmer)
        altstr (str/join alt)
        replacer #(str/replace % kmerstr altstr)]
    (reduce #(update-in %1 [%2 :seq] replacer)
            reads
            indices)))

(defn polish
  [k min-freq reads]
  (let [readseqs (map :seq reads)]
    (->> (build-spectrum k readseqs)
         (find-alternatives min-freq)
         (reduce fix-reads reads))))

(comment
  (polish 3 2 [{:seq "AGGTC"} {:seq "GTCTTGA"} {:seq "AAGGCTGTC"}])
  (def spectrum (build-spectrum 3 ["AGGTC" "GTCTTGA" "AAGGCTGTC"]))
  (def cs (candidates 2 spectrum))
  (apply max-key :freq cs)
  (find-alternatives 2 spectrum))
