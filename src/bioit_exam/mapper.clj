(ns bioit-exam.mapper)

(defn read-fasta [file-name])
  ; (with-open [r (BufferedReader. (FileReader.  file-name))]
  ; (.read )))

; (defn map-from-file [ref-file sample-file]
;   (let [reference (read-fasta ref-file)
;         reads (read-fasta sample-file)])
  ; (read-fasta ref-file)
  ; (read-fasta sample-file)

(defn get-seed [_read seed-len]
  (take seed-len _read))

(defn get-kmer-positions [_ref kmer]
  (let [seed-len (count kmer)])
  [1 2 3])

(defn pad-with [fill-in xs]
  (concat xs (repeat fill-in)))

(defn count-mismatches [ref' read' offset]
  (->> ref'
       (drop offset)
       (take (count read'))
       (pad-with nil)
       (map = read')
       (filter false?)
       (count)))

(count-mismatches "agtc" "gtac" 1)

(defn with-valid-positions [kmer-size max-mismatches _ref _read]
  (let [mismatches-tolerated?
        (fn [pos] (>= max-mismatches (count-mismatches _ref _read pos)))]
    (->> _read
         (get-seed kmer-size)
         (get-kmer-positions _ref)
         (filter mismatches-tolerated?)
         (reduce conj [])
         (conj [_read]))))

(defn map-reads [kmer-size max-mismatches _ref reads]
  (let [with-valid-positions (partial with-valid-positions kmer-size max-mismatches _ref)]
    (->> reads
         (map with-valid-positions)
         (into {}))))

(map-reads 23 2 "agtc" ["ag" "gt" "tc"])

