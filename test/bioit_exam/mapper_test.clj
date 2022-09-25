(ns bioit-exam.mapper-test
  (:require [clojure.test :refer [deftest testing is]]
            [bioit-exam.mapper :refer :all]))

(deftest map-reads-test
  (testing "should map single read at right position"
    (let [ref' "CTGACT"
          reads ["TGAC"]
          mappings (map-reads 2 2 ref' reads)]
      (is (= mappings {1 ["TGAC"]}))))

  (testing "should map all positions where kmers match"
    (let [ref' "CTGACTG"
          reads ["CTG"]
          mappings (map-reads 3 2 ref' reads)]
      (is (= mappings {0 ["CTG"], 4 ["CTG"]}))))

  (testing "should map reads when they have exactly the maximum of mismatches"
    (let [ref' "CTGACTG"
          reads ["TGAXTX"]
          mappings (map-reads 3 2 ref' reads)]
      (is (= mappings {1 ["TGAXTX"]}))))

  (testing "should not map reads when they exceed the maximum of mismatches"
    (let [ref' "CTGACTG"
          reads ["TGAXTX"]
          mappings (map-reads 3 3 ref' reads)]
      (is (= mappings {1 ["TGAXTX"]})))))
