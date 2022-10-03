(ns bioit-exam.spec
  (:require [clojure.spec.alpha :as s]))

(s/def ::base #(re-matches #"[ACGT]+" (str %)))

(s/def ::read (s/and string? ::base))

(s/def ::name string?)
(s/def ::seq (s/and vector?
                    (s/coll-of (s/and char? ::base))))
(s/def ::reference (s/keys :req-un [::name ::seq]))

(s/def ::mapped-reads
  (s/map-of nat-int? (s/coll-of ::read)))

(s/def ::mapping
  (s/keys :req-un [::reference ::mapped-reads]))
