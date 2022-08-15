(ns bioit-exam.io-test
  (:require [clojure.test :refer [deftest testing is]]
            [bioit-exam.io :refer :all]))

(deftest read-fasta-test
  (testing "should parse fasta lines properly"
    (let [reads (parse-fasta [">Read_1" "AGTCAGTC" "TCAGTCAG"
                              ">Read_2" "CCC" "AGAGAGA"])]
      (is (= 2 (count reads))))))

