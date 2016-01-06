(ns pressspan.test
  (:use [clojure.test])
  (:require [pressspan.graph :as pg]
            [pressspan.core :as pc]
            [pressspan.visualise :as pv]
            [pressspan.statistics :as ps]
            [pressspan.saminput :as psi]
            [clojure.data.int-map :as im]))



(deftest graph-clj-test
  (testing "graph.clj"
    (testing "structure"
      (let [A {:chr "1" :p5 1 :p3 177 :dir :plus}
            B {:chr "2" :p5 199 :p3 521 :dir :plus}
            C {:chr "2" :p5 1 :p3 127 :dir :plus}
            genome (ref {:frags (im/int-map), :nf 0, :links (im/int-map), :nl 0})]
        (pg/create-chr genome "1" 1000)
        (pg/create-chr genome "2" 1000)
        (pg/create-el genome {:chr "1" :p5 200 :p3 300 :dir :plus})
        (is (map? (pg/create-el genome A)))
        (is (map? (pg/create-el genome B)))
        (is (map? (pg/register-fragment genome (assoc C :prev B))))
        (is (map? (pg/get-el genome A)))
        (is (map? (pg/get-el genome B)))
        (is (nil? (pg/get-el genome (assoc B :dir :minus))))
        (is (nil? (pg/get-el genome (assoc A :chr "negative-test"))))
        (pg/link-frags genome (:id (pg/get-el genome A)) (:id (pg/get-el genome B)))
        (pg/link-frags genome (:id (pg/get-el genome A)) (:id (pg/get-el genome C)))))

    (testing ".sam file processing"
      (let [filename "test/data/5_out.sam"
            funs {:head? psi/header-line?
                  :header-parser psi/parse-header-line
                  :data-parser psi/make-frag
                  :add-all [pg/remember-multistrand pg/remember-circular pg/register-fragment]}
            genome (pg/create-genome filename funs)]
        (is (map? genome))
        (is (= 14 (:nf genome)))
        (is (= 10 (:nl genome)))
        (println (clojure.string/join \newline (pg/get-subgraph genome (first (:multis genome)) 1)))
        (is (= 4 (count (pg/get-subgraph genome (first (:multis genome))))))
        (is (= 5 (count (pg/all-subgraphs genome (:multis genome)))))))))

(deftest graph-test
  (testing "visualise.clj"
    (let [filename "test/data/5_out.sam"
          funs {:head? psi/header-line?
                :header-parser psi/parse-header-line
                :data-parser psi/make-frag
                :add-all [pg/remember-multistrand pg/register-fragment]}
          genome (pg/create-genome filename funs)
          multi (first (:multis genome))
          subgraph (pg/get-subgraph genome multi)]
      (is (true? ((pv/drop-filter 2) [{:up [{:depth 2} {:depth 5}]} {:up [{:depth 3}]} {:down nil}])))
      (is (false? ((pv/drop-filter 3) [{:up [{:depth 2} {:depth 5}]} {:up [{:depth 3}]} {:down nil}])))
      (pv/write-files genome :multis "test/output" 1)
      (println (pv/graph->dot genome subgraph "testgraph"))
      (println (pv/graph->log genome subgraph "testgraph")))))


(deftest stat-test
  (testing "statistics.clj"
    (let [filename "test/data/5_out.sam"
          funs {:head? psi/header-line?
                :header-parser psi/parse-header-line
                :data-parser psi/make-frag
                :add-all [pg/remember-multistrand pg/register-fragment]}
          genome (pg/create-genome filename funs)]
      (ps/empty-stats genome 100)
      (ps/get-pairs genome ps/is-multi? (:multis genome))
      (ps/write-stat-file "test/output/multi.csv" :multis 1000 genome))))
