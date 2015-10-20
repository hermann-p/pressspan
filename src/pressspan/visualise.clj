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
      (trampoline find-here (set graph) #{} #{} vert-map {} 0 ))))


(defn insert-dummies [graph]
  (let [did -1]
    ))


(defn- fix-links
	"Changes el->link-id->link->el-id structure to el->el-id"
	[el root]
	(let [up-links (map #(get-in root [:links %]) (:up el))
;			  up-ids (map :up up-links)
			  dn-links (map #(get-in root [:links %]) (:dn el))]
	;		  dn-ids (map :dn dn-links)]
		(-> el
	  	(assoc :up up-links)
	  	(assoc :down dn-links))))


(defn- load-and-fix
  "Transforms a list of element-ids to a map of element-id -> element with 
  link-ids changed to linked-element-ids"
  [graph root]
  (let [graph (map #(get-in root [:frags %]) graph)]
  (zipmap (map :id graph) (map #(fix-links % root) graph))))


(defn- hex [n]
  {:pre [(integer? n) (>= n 0) (<= n 255)]}
  (if (> n 16)
    (format "%x" n)
    (format "0%x" n)))

(defn- chr-color [n]
  (nth (cycle (for [r (map (partial * 85) (range 4))
             g [255 0 128]
             b (map #(- 255 (* 85 %)) (range 4))]
         (str "#" (hex r) (hex g) (hex b))))
       (- (* 3 n) 1)))


(defn graph->dot [root graph id-string]
	(let [g (vals (load-and-fix graph root))
				put-el
				(fn [el]
          (cons 
            (str "  " (:id el) " [shape=\"rectangle\",border=0,style=\"filled\",height=0.25,"
                 "label=\"" (:chr el) ":" (:p5 el) "-" (:p3 el) "\","
                 "color=\"" (chr-color (Integer. (:chr el))) "\"];\n")
            (for [link (:up el)] (str "  " (:id el) "->" (:up link) " [label=\"" (:depth link) "\"];\n"))))
	      walk
	      (fn walk [els]
	      	(if (seq els)
	         (conj (walk (rest els)) (clojure.string/join (put-el (first els))))))]
	  (str "digraph " id-string " {" \newline 
         "  rankdir=\"LR\";" \newline
	       (clojure.string/join (doall (concat (trampoline walk g))))
	       "}")))


(defn graph->log [root graph id-string]
  (let [nodes (map #(get-in root [:frags %]) graph)]
    (format "%s\t|%s|\t%s\n"
    	id-string 
	    (clojure.string/join ","
	      (for [el nodes] 
	        (str (if (= :minus (:dir el)) "-") (:chr el) ":" (:p5 el) "-" (:p3 el))))
      (clojure.string/join ","
	      (sort (set (map :chr nodes)))))))


(defn write-files [root type basedir]
  (let [graphs (pressspan.graph/all-subgraphs root (type root))
  		  typestr (name type)]
    (println "Writing" (count graphs) typestr "graph files")
    (pressspan.io/make-dir (str basedir "/" typestr))
    (doall
      (map-indexed
        #(spit
            (str basedir "/" typestr "/" typestr "_" %1 ".dot")
            (graph->dot root %2 (str typestr "_" %1)))
        graphs)))
  root)
                    


(deftest graph-test
  (let [filename "test/data/5_out.sam"
        funs {:head? pressspan.saminput/header-line?
              :header-parser pressspan.saminput/parse-header-line
              :data-parser pressspan.saminput/make-frag
              :add-all [pressspan.graph/remember-multistrand pressspan.graph/register-fragment]}
        genome (pressspan.graph/create-genome filename funs)
        multi (first (:multis genome))
        subgraph (pressspan.graph/get-subgraph genome multi)]
    (is (map? (assign-layers genome subgraph)))
    (is (map? (fix-links (get-in genome [:frags (first subgraph)]) genome)))
    (is (= 4 (count (load-and-fix subgraph genome))))
    (write-files genome :multis "test/output")
    (println (graph->dot genome subgraph "testgraph"))
    (println (graph->log genome subgraph "testgraph"))))