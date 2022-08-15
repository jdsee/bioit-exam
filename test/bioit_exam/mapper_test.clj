(ns bioit-exam.mapper-test
  (:require [clojure.test :refer [deftest testing is]]
            [bioit-exam.mapper :refer :all]))

(deftest map-reads-test
  (testing "should map single read at right position"
    (let [ref' "CTGACT"
          reads ["TGAC"]
          mappings (map-reads 2 2 ref' reads)]
      (is (= mappings {"TGAC" [1]}))))

  (testing "should map all positions where kmers match"
    (let [ref' "CTGACTG"
          reads ["CTG"]
          mappings (map-reads 3 2 ref' reads)]
      (is (= mappings {"CTG" [0 4]}))))

  (testing "should map reads when they have exactly the maximum of mismatches"
    (let [ref' "CTGACTG"
          reads ["TGAXTX"]
          mappings (map-reads 3 2 ref' reads)]
      (is (= mappings {"TGAXTX" [1]}))))

  (testing "should not map reads when they exceed the maximum of mismatches"
    (let [ref' "CTGACTG"
          reads ["TGAXTX"]
          mappings (map-reads 3 3 ref' reads)]
      (is (= mappings {"TGAXTX" [1]})))))

