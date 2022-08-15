(ns bioit-exam.reader-test
  (:require [clojure.test :refer [deftest testing is]]
            [bioit-exam.reader :refer :all]))

(deftest read-fasta-test
  (testing "should parse fasta lines properly"
    (let [reads (parse-fasta [">Read_1" "AGTCAGTC" "TCAGTCAG"
                              ">Read_2" "CCC" "AGAGAGA"])]
      (is (= [{:header ">Read_1" :lines ["AGTCAGTC" "TCAGTCAG"]}
              {:header ">Read_2" :lines ["CCC" "AGAGAGA"]}]
             reads)))))

