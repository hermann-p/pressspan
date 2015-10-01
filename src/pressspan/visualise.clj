(ns pressspan.visualise
  (:use [clojure.test]
        [clojure.set :only [difference union]])
  (:require [pressspan.graph :refer [all-subgraphs]]))


(defn assign-layers
  "returns graph element with assigned layers.
   root  - map create by pressspan.graph/create-genome
   graph - seq of element indices"
  [root graph]
  (letfn
    [(up-known? [up-known verts el-id]
                (let [el (get verts el-id)
                      links (map #(get-in root [:links %]) (:up el))
                      this-up (set (map :up links))]
                  (empty? (difference this-up up-known))))    ; el has no unknown upstream links
     
     (find-next [free verts up-els]
                (first (filter (partial up-known? up-els verts) free)))
     
     (find-here [free up-els els-here vert-map layered-graph layer]
                (if (seq free)
                  (if-let [v (find-next free vert-map up-els)]
                               
                    (find-here                            ; map found el to current layer
                      (disj free v)                            ; remove it from available els
                      up-els                                   ; leave upper layer elements alone
                      (conj els-here v)                        ; store el in els of this layer
                      (assoc-in vert-map [v :layer] layer)     ; assign layer number to el
                      (update-in layered-graph [layer] conj v) ; add element to layered graph
                      layer)                                   ; keep layer number
                    (find-here                            ; step up one layer
                      free                                     ; keep free elements
                      (union up-els els-here)                  ; save elements in "known from upper layer"
                      #{}                                      ; empty set of els in new layer
                      vert-map                                 ; keep vertices
                      layered-graph                            ; keep layered graph
                      (inc layer)))                            ; increase layer number
                  layered-graph))]
    (let [vertices (mapv #(get-in root [:frags %]) graph)
          vert-map (zipmap graph vertices)]
      (println "vert-map:")
      (trampoline find-here (set graph) #{} #{} vert-map {} 0 ))))


(defn insert-dummies [graph]
  (let [did -1]
    ))


(deftest graph-test
  (let [filename "test/data/5_out.sam"
        funs {:head? pressspan.saminput/header-line?
              :header-parser pressspan.saminput/parse-header-line
              :data-parser pressspan.saminput/make-frag
              :add-all [pressspan.graph/remember-multistrand pressspan.graph/register-fragment]}
        genome (pressspan.graph/create-genome filename funs)
        multi (first (:multis genome))
        subgraph (pressspan.graph/get-subgraph genome multi)]
    (is (map? (assign-layers genome subgraph)))))