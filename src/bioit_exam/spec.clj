(ns bioit-exam.spec
  (:require [clojure.spec.alpha :as s]))

(s/def ::read
  (s/and string? #(re-matches #"[ACGT]+" %)))

(s/def ::reference ::read)

(s/def ::mapped-reads
  (s/map-of nat-int? (s/coll-of ::read)))

(s/def ::mapping
  (s/and (s/keys :req [::reference ::mapped-reads])))

