(ns pressspan.first-tests
  (:require pressspan.types
            pressspan.core)
  (:import [pressspan.types Genome Chromosome Fragment]))

(use 'clojure.test)

(deftest datatypes
  (is (=
       (-> (add-chromosome (Genome. '()) "A" 1000) :chrs (nth 0) :id-str)
       "A"))
  
  (is (=
       (-> (add-fragment (Chromosome. 5 "x" 1000 '()) 100 110) :fragments (nth 0) :pos3)
       110)))
