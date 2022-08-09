(ns bioit-exam.core
  (:require [clojure.spec.alpha :as s]
            [clojure.string :as str]))

(s/def ::base (s/and string? ))
(s/def ::bases (s/coll-of ::base))

