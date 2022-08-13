(ns bioit-exam.polisher-test
  (:require [clojure.test :refer [deftest testing is]]
            [bioit-exam.polisher :refer :all]))

(deftest polish-kmer-test
  (testing "polishing replaces char in the middle of a kmer"
    (is (= (polish "AGTC" "X" 2) [\A \G \X \C])))

  (testing "polishing replaces char at the front of a kmer"
    (is (= (polish "AGTC" "X" 0) [\X \G \T \C])))

  (testing "polishing replaces char at the end of a kmer"
    (is (= (polish "AGTC" "X" 3) [\A \G \T \X]))))

(deftest build-spectrum-test
  (testing "building spectrum with correct frequencies"
    (is (= (build-spectrum 2 ["ATAGTC" "GTCATC"])
           {[\T \C] {:reads ["ATAGTC" "GTCATC"], :freq 3},
            [\A \T] {:reads ["ATAGTC" "GTCATC"], :freq 2},
            [\G \T] {:reads ["ATAGTC" "GTCATC"], :freq 2},
            [\T \A] {:reads ["ATAGTC"], :freq 1},
            [\A \G] {:reads ["ATAGTC"], :freq 1},
            [\C \A] {:reads ["GTCATC"], :freq 1}}))))

(deftest find-best-candidate-test
  (testing "choosing candidates is based on the highest frequencies in spectrum"
    (let [chosen-candidates (choose-candidates
                             3
                             {"AGT" {:freq 2}, "GGT" {:freq 3},"AGC" {:freq 5}}
                             [{:kmer "AGT" :replacement "GGT"}
                              {:kmer "AGT" :replacement "AGC"}])
          chosen-replacements (map :replacement chosen-candidates)]
      (is (.contains chosen-replacements "AGC"))))

  (testing "choosing candidates does not crash when replacement is missing in spectrum"
    (let [chosen-candidates (choose-candidates
                             3
                             {"AGT" {:freq 2}, "AGC" {:freq 5}}
                             [{:kmer "AGT" :replacement "AGC"}
                              {:kmer "AGT" :replacement "XXX"}])]
      (is (some? (vec chosen-candidates)))))

  (testing "choosing candidates returns one replacement for each kmer at max"
    (let [chosen-candidates (choose-candidates
                             3
                             {"AGT" {:freq 2}, "GGT" {:freq 3}, "GCC" {:freq 7}, "AGC" {:freq 5}}
                             [{:kmer "AGT" :replacement "GGT"}
                              {:kmer "AGT" :replacement "AGC"}
                              {:kmer "GTC" :replacement "GCC"}])
          chosen-replacements (map :replacement chosen-candidates)]
      (is (= chosen-replacements ["AGC" "GCC"]))))

  (testing "choosing candidates return not more than one replacement per kmer"
    (let [chosen-candidates (choose-candidates
                             3
                             {"AGT" {:freq 2}, "GGT" {:freq 3}, "AGC" {:freq 5}}
                             [{:kmer "AGT" :replacement "GGT"}
                              {:kmer "AGT" :replacement "AGC"}])]
      (is (= 1 (count chosen-candidates)))))

  (testing "choosing candidates does not return replacements with frequency lower than the minimum"
    (let [chosen-candidates (choose-candidates
                             3
                             {"AGT" {:freq 2}, "GGT" {:freq 3}, "AGC" {:freq 5}}
                             [{:kmer "AGT" :replacement "GGT"}
                              {:kmer "AGT" :replacement "AGC"}])]
      (is (= 1 (count chosen-candidates))))))

